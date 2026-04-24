package com.mcgui.javagui;

import javafx.fxml.FXML;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.*;
import javafx.scene.control.ScrollPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.application.Platform;
import java.util.List;
import java.util.ArrayList;
import javafx.scene.text.Text;
import javafx.scene.chart.NumberAxis;
import javafx.scene.layout.Pane;
import javafx.scene.layout.AnchorPane;
import javafx.scene.shape.Line;
import javafx.scene.paint.Color;

public class GUIController {
    @FXML
    private Button startSimButton;
    @FXML
    private Button adLayerButton;
    @FXML
    private ScrollPane layerScroll;
    @FXML
    private VBox layerContainer;
    @FXML private LineChart<Number, Number> trajectoryChart;
    @FXML private TextField photonsField;
    @FXML private ProgressBar progressBar;
    @FXML private Label progressLabel;
    @FXML private NumberAxis xAxis;
    @FXML private NumberAxis yAxis;
    private Pane boundaryPane;  // контейнер для линий-границ

    private MCSimulationJNI jni = new MCSimulationJNI();
    private List<LayerParams> lastLayers = null;

    @FXML
    public void initialize() {
        adLayerButton.setOnAction(event->ShowAddLayerDialoge());
        startSimButton.setOnAction(event -> onRunSimulation());

        trajectoryChart.lookup(".chart-line-symbol").setStyle(
                "-fx-background-radius: 0.1px; -fx-padding: 0.1px; -fx-background-color: rgba(100,100,100,0.3);"
        );

        progressBar.setProgress(0.0);
        progressLabel.setText("Готов к запуску");

        boundaryPane = new Pane();
        boundaryPane.setMouseTransparent(true);  // чтобы не мешать взаимодействию с графиком
        AnchorPane parent = (AnchorPane) trajectoryChart.getParent();
        parent.getChildren().add(boundaryPane);
// Привязываем к тем же якорям, что и график
        AnchorPane.setTopAnchor(boundaryPane, AnchorPane.getTopAnchor(trajectoryChart));
        AnchorPane.setBottomAnchor(boundaryPane, AnchorPane.getBottomAnchor(trajectoryChart));
        AnchorPane.setLeftAnchor(boundaryPane, AnchorPane.getLeftAnchor(trajectoryChart));
        AnchorPane.setRightAnchor(boundaryPane, AnchorPane.getRightAnchor(trajectoryChart));

// Обновляем границы при изменении размеров графика
        trajectoryChart.widthProperty().addListener((obs, old, val) -> updateBoundaries());
        trajectoryChart.heightProperty().addListener((obs, old, val) -> updateBoundaries());
        yAxis.lowerBoundProperty().addListener((obs, old, val) -> updateBoundaries());
        yAxis.upperBoundProperty().addListener((obs, old, val) -> updateBoundaries());
    }

    private void ShowAddLayerDialoge(){
        Dialog<LayerParams> dialog = new Dialog<>();
        dialog.setTitle("Добавлениие слоя");
        dialog.setHeaderText("Введите параметры слоя");

        ButtonType okButtonType = new ButtonType("Добавить", ButtonBar.ButtonData.OK_DONE);
        dialog.getDialogPane().getButtonTypes().addAll(okButtonType, ButtonType.CANCEL);

        TextField muAField = new TextField();
        muAField.setPromptText("μa");
        TextField muSField=new TextField();
        muSField.setPromptText("μs");
        TextField gField = new TextField();
        gField.setPromptText("g");
        TextField nField = new TextField();
        nField.setPromptText("n");
        TextField z_botField = new TextField();
        z_botField.setPromptText("z_bot");
        TextField z_topField = new TextField();
        z_topField.setPromptText("z_top");

        GridPane grid = new GridPane();
        grid.setHgap(10);
        grid.setVgap(10);
        grid.setPadding(new Insets(20, 150, 10, 10));
        grid.add(new Label("μa"),0,0);
        grid.add(muAField,1,0);
        grid.add(new Label("μs"),0,1);
        grid.add(muSField,1,1);
        grid.add(new Label("g"),0,2);
        grid.add(gField,1,2);
        grid.add(new Label("n"),0,3);
        grid.add(nField,1,3);
        grid.add(new Label("z_bot"),0,4);
        grid.add(z_botField,1,4);
        grid.add(new Label("z_top"),0,5);
        grid.add(z_topField,1,5);
        dialog.getDialogPane().setContent(grid);

        dialog.setResultConverter(dialogButton -> {
            if (dialogButton == okButtonType) {
                try {
                    double mua = Double.parseDouble(muAField.getText());
                    double mus = Double.parseDouble(muSField.getText());
                    double g = Double.parseDouble(gField.getText());
                    double n = Double.parseDouble(nField.getText());
                    double z_top = Double.parseDouble(z_topField.getText());
                    double z_bot = Double.parseDouble(z_botField.getText());
                    return new LayerParams(mua, mus, g, n, z_top, z_bot);
                } catch (NumberFormatException e) {
                    return null;
                }
            }
            return null;
        });
        dialog.showAndWait().ifPresent(params -> addLayerToUI(params));
    }
    private void addLayerToUI(LayerParams params) {
        // Создаём панель для отображения слоя
        HBox layerBox = new HBox(15);
        layerBox.setAlignment(Pos.CENTER_LEFT);
        layerBox.setStyle("-fx-border-color: lightgray; -fx-padding: 5; -fx-background-color: #f9f9f9;");

        layerBox.setUserData(params);
        Label infoLabel = new Label(String.format("μa=%.4f  μs=%.2f  g=%.2f  n=%.2f  z_top=%.2f  z_bot=%.2f",
                params.mua, params.mus, params.g, params.n, params.z_top, params.z_bot));
        infoLabel.setPrefWidth(400);

        Button removeButton = new Button("Удалить");
        removeButton.setOnAction(e -> layerContainer.getChildren().remove(layerBox));

        layerBox.getChildren().addAll(infoLabel, removeButton);
        layerContainer.getChildren().add(layerBox);
    }
    @FXML
    private void onRunSimulation() {
        startSimButton.setDisable(true);
        progressBar.setProgress(0.0);
        progressLabel.setText("Запуск...");

        // Читаем число фотонов
        long numPhotons;
        try {
            numPhotons = Long.parseLong(photonsField.getText());
            if (numPhotons <= 0) throw new NumberFormatException();
        } catch (NumberFormatException e) {
            progressLabel.setText("Ошибка: введите положительное число фотонов");
            startSimButton.setDisable(false);
            return;
        }
        List<LayerParams> layers = new ArrayList<>();
        for (Node node : layerContainer.getChildren()) {
            if (node instanceof HBox box) {
                // Параметры уже сохранены в Label, но лучше хранить их в поле каждого HBox.
                // Упростим: будем хранить LayerParams как свойство HBox.
                Object tag = box.getUserData();
                if (tag instanceof LayerParams lp) layers.add(lp);
            }
        }
        if (layers.isEmpty()) {
            progressLabel.setText("Ошибка: добавьте хотя бы один слой");
            startSimButton.setDisable(false);
            return;
        }
        this.lastLayers = new ArrayList<>(layers);
        // Подготовить плоский массив double
        double[] flat = new double[layers.size() * 6];
        for (int i = 0; i < layers.size(); i++) {
            LayerParams lp = layers.get(i);
            flat[i*6]   = lp.mua;
            flat[i*6+1] = lp.mus;
            flat[i*6+2] = lp.g;
            flat[i*6+3] = lp.n;
            flat[i*6+4] = lp.z_top;
            flat[i*6+5] = lp.z_bot;
        }
        // Запустить в отдельном потоке, чтобы не блокировать UI
        new Thread(() -> {
            Platform.runLater(() -> {
                progressBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
                progressLabel.setText("Выполняется симуляция...");
            });

            // Вызов JNI (блокирующий)
            Object[][] trajs = jni.runSimulation(flat, numPhotons, 0,0,0, 0,0,1, 1.0);

            Platform.runLater(() -> {
                progressBar.setProgress(1.0);
                progressLabel.setText("Симуляция завершена");
                drawTrajectories(trajs);
                startSimButton.setDisable(false);
            });
        }).start();
    }
    private void drawTrajectories(Object[][] trajectories) {
        // Очищаем график
        trajectoryChart.getData().clear();

        for (Object trajObj : trajectories) {
            double[] points = (double[]) trajObj;
            if (points.length < 6) continue; // нужно минимум 2 точки

            // Новая серия для каждой траектории
            XYChart.Series<Number, Number> series = new XYChart.Series<>();
            //series.setName("Фотон " + (trajectoryChart.getData().size() + 1));

            // Добавляем точки в порядке следования (они будут соединены линиями)
            for (int i = 0; i < points.length / 3; i++) {
                double x = points[i * 3];
                double z = points[i * 3 + 2]; // вид XZ
                series.getData().add(new XYChart.Data<>(x, z));
            }

            trajectoryChart.getData().add(series);
            drawLayerBoundaries(this.lastLayers);
        }
    }
    private void drawLayerBoundaries(List<LayerParams> layers) {
        boundaryPane.getChildren().clear();
        for (LayerParams lp : layers) {
            addBoundaryLine(lp.z_top);
            addBoundaryLine(lp.z_bot);
        }
    }

    private void addBoundaryLine(double z) {
        double y = yAxis.getDisplayPosition(z);
        if (Double.isNaN(y) || y < 0 || y > trajectoryChart.getHeight()) return;
        Line line = new Line(0, y, trajectoryChart.getWidth(), y);
        line.setStroke(Color.LIGHTGRAY);
        line.setStrokeWidth(0.8);
        line.getStrokeDashArray().addAll(5.0, 5.0);  // пунктир
        boundaryPane.getChildren().add(line);
    }

    private void updateBoundaries() {
        // вызывается при изменении размеров или диапазона осей
        if (lastLayers != null) drawLayerBoundaries(lastLayers);
    }
}

