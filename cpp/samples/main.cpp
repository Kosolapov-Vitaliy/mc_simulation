#include <iostream>
#include <rund_num_generate.h>
#include <clocale>
int main() {
    std::setlocale(LC_ALL, "Russian");
    double av_len = 0.5;
    RNGenerate generator;
    for (int i = 0; i < 5; i++) {
        double temp = generator.LengthGenerate(av_len);
        std::cout <<"Шаг:" << temp<<"\n";
    }
    std::cout << "Закончите программу";
    std::cin >> av_len;
}
