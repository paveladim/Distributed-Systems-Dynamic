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