#include <iostream>
#include <vector>


// 1. Создать папку сборки
// mkdir build
// cd build
// 2. Сгенерировать проект
// cmake ..
// 3. Собрать
// make

using state_t = std::tuple<size_t, size_t, size_t>;
state_t foo(int a, int b, int c){
    return std::make_tuple(a+b, b, c-1);
}
int main() {
    int a, b, c;
    a = 2;
    b = 3;
    c = 4;
    std::tie(a,b,c) = foo(a, b, c);
    std::cout << a << " " << b << " " << c << std::endl;

    int x = -5;
    int y = -3;
    std::cout << x%y << std::endl;
    return 0;
}