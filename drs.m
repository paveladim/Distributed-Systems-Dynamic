function drs
    xlist = 7.5 : 0.001 : 12;
    
    f = @(x)(cos(x) * sinh(x) - sin(x) * cosh(x));
    
    flist = [];
    for i = 1:length(xlist)
        flist = [flist, f(xlist(i))];
    end

    x0 = [3 4];
    b1 = fzero(f, x0);
    x0 = [7 8];
    b2 = fzero(f, x0);
    x0 = [10 11];
    b3 = fzero(f, x0);

    C11 = 1;
    C12 = sin(b1) / sinh(b1);

    C21 = 1;
    C22 = sin(b2) / sinh(b2);

    C31 = 1;
    C32 = sin(b3) / sinh(b3);

    ksilist = 0 : 0.001 : 1;

    v1 = @(x)(C11 * sin(b1 * x) + C12 * sinh(b1 * x));
    v2 = @(x)(C21 * sin(b2 * x) + C22 * sinh(b2 * x));
    v3 = @(x)(C31 * sin(b3 * x) + C32 * sinh(b3 * x));

    v12 = @(x)((C11 * sin(b1 * x) + C12 * sinh(b1 * x))^2);
    norm1 = sqrt(integral(v12, 0, 1, "ArrayValued", true));
    v22 = @(x)((C21 * sin(b2 * x) + C22 * sinh(b2 * x))^2);
    norm2 = sqrt(integral(v22, 0, 1, "ArrayValued", true));
    v32 = @(x)((C31 * sin(b3 * x) + C32 * sinh(b3 * x))^2);
    norm3 = sqrt(integral(v32, 0, 1, "ArrayValued", true));

    v1l = [];
    v2l = [];
    v3l = [];
    for i = 1:length(ksilist)
        v1l = [v1l, v1(ksilist(i)) / norm1];
        v2l = [v2l, v2(ksilist(i)) / norm2];
        v3l = [v3l, v3(ksilist(i)) / norm3];
    end

    x0 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    x0 = fsolve(@syste, x0);

    fhandle = figure;
    subplot(3, 1, 1)
        plot(ksilist, v1l, 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 2)
        plot(ksilist, v2l, 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 3)
        plot(ksilist, v3l, 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');

end

function F = syste(X)
    F(1) = 6 * X(2) + 12 * X(3);
    F(2) = 6 * X(5) + 12 * X(6) + 20 * X(7);
    F(3) = 6 * X(9) + 12 * X(10) + 20 * X(11) + 30 * X(12);
    F(4) = 6 * X(2) + 24 * X(3);
    F(5) = 6 * X(5) + 24 * X(6) + 60 * X(7);
    F(6) = 6 * X(9) + 24 * X(10) + 60 * X(11) + 120 * X(12);

    F(7) = X(1) * X(1) / 3 + 2 * X(1) * X(2) / 5 + 2 * X(1) * X(3) / 6;
    F(7) = F(7) + X(2) * X(2) / 7 + 2 * X(2) * X(3) / 8 + X(3) * X(3) / 9 - 1;

    F(8) = X(4) * X(4) / 3 + 2 * X(4) * X(5) / 5 + 2 * X(4) * X(6) / 6 + 2 * X(4) * X(7) / 7 + X(5) * X(5) / 7;
    F(8) = F(8) + 2 * X(5) * X(6) / 8 + 2 * X(5) * X(7) / 9 + X(6) * X(6) / 9 + 2 * X(6) * X(7) / 10 + X(7) * X(7) / 11 - 1;

    F(9) = X(8) * X(8) / 3 + 2 * X(8) * X(9) / 5 + 2 * X(8) * X(10) / 6 + 2 * X(8) * X(11) / 7 + 2 * X(8) * X(12) / 8;
    F(9) = F(9) + X(9) * X(9) / 7 + 2 * X(9) * X(10) / 8 + 2 * X(9) * X(11) / 9 + 2 * X(9) * X(12) / 10;
    F(9) = F(9) + X(10) * X(10) / 9 + 2 * X(10) * X(11) / 10 + 2 * X(10) * X(12) / 11;
    F(9) = F(9) + X(11) * X(11) / 11 + 2 * X(11) * X(12) / 12 + X(12) * X(12) / 13 - 1;

    F(10) = X(1) * X(4) / 3 + X(1) * X(5) / 5 + X(1) * X(6) / 6 + X(1) * X(7) / 7 + X(2) * X(4) / 5 + X(2) * X(5) / 7;
    F(10) = F(10) + X(2) * X(6) / 8 + X(2) * X(7) / 9 + X(3) * X(4) / 6 + X(3) * X(5) / 8 + X(3) * X(6) / 9 + X(3) * X(7) / 10;

    F(11) = X(1) * X(8) / 3 + X(1) * X(9) / 5 + X(1) * X(10) / 6 + X(1) * X(11) / 7 + X(1) * X(12) / 8 + X(2) * X(8) / 5 + X(2) * X(9) / 7;
    F(11) = F(11) + X(2) * X(10) / 8 + X(2) * X(11) / 9 + X(2) * X(12) / 10 + X(3) * X(8) / 6 + X(3) * X(9) / 8 + X(3) * X(10) / 9;
    F(11) = F(11) + X(3) * X(11) / 10 + X(3) * X(12) / 11;

    F(12) = X(4) * X(8) / 3 + X(4) * X(9) / 5 + X(4) * X(10) / 6 + X(4) * X(11) / 7 + X(4) * X(12) / 8;
    F(12) = F(12) + X(5) * X(8) / 5 + X(5) * X(9) / 7 + X(5) * X(10) / 8 + X(5) * X(11) / 9 + X(5) * X(12) / 10;
    F(12) = F(12) + X(6) * X(8) / 6 + X(6) * X(9) / 8 + X(6) * X(10) / 9 + X(6) * X(11) / 10 + X(6) * X(12) / 11;
    F(12) = F(12) + X(7) * X(8) / 7 + X(7) * X(9) / 9 + X(7) * X(10) / 10 + X(7) * X(11) / 11 + X(7) * X(12) / 12;
end