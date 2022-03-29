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

    f11 = @(x)(b1^2 * (-C11 * sin(b1 * x) + C12 * sinh(b1 * x)) * (C11 * sin(b1 * x) + C12 * sinh(b1 * x)) / (norm1 * norm1));
    f12 = @(x)(b2^2 * (-C21 * sin(b2 * x) + C22 * sinh(b2 * x)) * (C11 * sin(b1 * x) + C12 * sinh(b1 * x)) / (norm1 * norm2));

    f21 = @(x)(b1^2 * (-C11 * sin(b1 * x) + C12 * sinh(b1 * x)) * (C21 * sin(b2 * x) + C22 * sinh(b2 * x)) / (norm1 * norm2));
    f22 = @(x)(b2^2 * (-C21 * sin(b2 * x) + C22 * sinh(b2 * x)) * (C21 * sin(b2 * x) + C22 * sinh(b2 * x)) / (norm2 * norm2));

    a11 = integral(f11, 0, 1, "ArrayValued", true);
    a21 = integral(f12, 0, 1, "ArrayValued", true);

    a12 = integral(f21, 0, 1, "ArrayValued", true);
    a22 = integral(f22, 0, 1, "ArrayValued", true);

    v1l = [];
    v2l = [];
    v3l = [];
    for i = 1:length(ksilist)
        v1l = [v1l, v1(ksilist(i)) / norm1];
        v2l = [v2l, v2(ksilist(i)) / norm2];
        v3l = [v3l, v3(ksilist(i)) / norm3];
    end

    x0 = ones(1, 12);
    x0 = fsolve(@syste, x0);

    D = -ones(1, 6);
    fun = @(x)syst(x, x0);
    D = fsolve(fun, D);

    p1l = [];
    p2l = [];
    p3l = [];
    for i = 1:length(ksilist)
        ksi = ksilist(i);
        p1l = [p1l, x0(4) * ksi + x0(5) * ksi^3 + x0(6) * ksi^4 + x0(7) * ksi^5];
        p2l = [p2l, x0(8) * ksi + x0(9) * ksi^3 + x0(10) * ksi^4 + x0(11) * ksi^5 + x0(12) * ksi^6];
        p3l = [p3l, D(1) * ksi + D(2) * ksi^3 + D(3) * ksi^4 + D(4) * ksi^5 + D(5) * ksi^6 + D(6) * ksi^7];
    end

    disp1 = 0;
    disp2 = 0;
    disp3 = 0;
    for i = 2:length(ksilist)
        disp1 = disp1 + (p1l(i) - v1l(i))^2;
        disp2 = disp2 + (p2l(i) - v2l(i))^2;
        disp3 = disp3 + (p3l(i) - v3l(i))^2;
    end

    disp1 = disp1 * 0.001;
    disp2 = disp2 * 0.001;
    disp3 = disp3 * 0.001;

    A0 = 1;
    A1 = @(b)(b * (a11 + a22) + b1^4 + b2^4);
    A2 = @(b)((b * a11 + b1^4) * (b * a22 + b2^4) - b^2 * a12 * a21);
    
    blst = 20.15 : 0.000001 : 20.2;
    unstablelst11 = [];
    unstablelst12 = [];
    unstablelst2 = [];

    b_crit = 0;
    for i = 1:length(blst)
        el = blst(i);
        unstablelst11 = [unstablelst11; -A2(el) + sqrt(A1(el) * A1(el) - 4 * A2(el))];
        unstablelst12 = [unstablelst12; -A2(el) - sqrt(A1(el) * A1(el) - 4 * A2(el))];
        if (unstablelst11(end) > 0) && (b_crit == 0)
            b_crit = el;
        end
        unstablelst2 = [unstablelst2; A1(el) * A1(el) - 4 * A2(el)];
    end

    b_crit;

    fhandle = figure;
    subplot(3, 2, 1)
        plot(ksilist, v1l, 'r', ksilist, p1l, 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 2, 2)
        plot(ksilist, v2l, 'r', ksilist, p2l, 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 2, 3)
        plot(ksilist, v3l, 'r', ksilist, p3l, 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 2, 4)
        plot(blst, unstablelst2, 'g', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 2, 5)
        plot(blst, unstablelst11, 'g', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 2, 6)
        plot(blst, unstablelst12, 'g', 'LineWidth', 2.0)
        grid on;
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('f(x)', 'FontSize', 12, 'FontWeight', 'bold');

end

function F = syst(X, coef)
    a = zeros(1, 5);
    b = zeros(1, 6);
    c = zeros(1, 7);
    d = zeros(1, 8);

    a(2) = coef(1);
    a(4) = coef(2);
    a(5) = coef(3);

    b(2) = coef(4);
    b(4) = coef(5);
    b(5) = coef(6);
    b(6) = coef(7);

    c(2) = coef(8);
    c(4) = coef(9);
    c(5) = coef(10);
    c(6) = coef(11);
    c(7) = coef(12);

    d(2) = X(1);
    d(4) = X(2);
    d(5) = X(3);
    d(6) = X(4);
    d(7) = X(5);
    d(8) = X(6);

    F = zeros(1, 6);

    for i = 2:7
        F(1) = F(1) + i * (i - 1) * d(i + 1);
    end

    for i = 3:7
        F(2) = F(2) + i * (i - 1) * (i - 2) * d(i + 1);
    end

    for i = 0:7
        for j = 0:7
            F(3) = F(3) + d(i + 1) * d(j + 1) / (i + j + 1);
        end
    end

    F(3) = F(3) - 1;

    for i = 0:7
        for j = 0:6
            F(4) = F(4) + d(i + 1) * c(j + 1) / (i + j + 1);
        end
    end

    for i = 0:7
        for j = 0:5
            F(5) = F(5) + d(i + 1) * b(j + 1) / (i + j + 1);
        end
    end

    for i = 0:7
        for j = 0:4
            F(6) = F(6) + d(i + 1) * a(j + 1) / (i + j + 1);
        end
    end
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