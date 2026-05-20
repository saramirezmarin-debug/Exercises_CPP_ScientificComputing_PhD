#include <iostream>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/color.h>

using Eigen::MatrixXd;

int main()
{
    // =========================
    // Eigen
    // =========================
    MatrixXd m(2,2);

    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);

    // =========================
    // fmt formatting
    // =========================
    int a = 16165;
    int b = 161634;

    std::string res;

    res = fmt::format(
        "The sum of {:^20} and {:-^20} is {}\n",
        a,
        b,
        a + b
    );

    fmt::print("{}\n", res);

    // =========================
    // Matrix printing
    // =========================
    fmt::print(
        fg(fmt::color::cyan),
        "Eigen matrix:\n"
    );

    std::cout << m << std::endl;

    // =========================
    // Colored output
    // =========================
    fmt::print(
        fg(fmt::color::red),
        "\nThis is red text\n"
    );

    fmt::print(
        bg(fmt::color::green),
        "This is green background\n"
    );

    fmt::print(
        fg(fmt::color::blue) | bg(fmt::color::yellow),
        "Blue text on yellow background\n"
    );

    return 0;
}