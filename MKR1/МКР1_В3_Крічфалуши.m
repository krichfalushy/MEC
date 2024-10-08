%  --- Задача #1 ---


A = [20, 2, -2, 4, -3.25;
     2, 15, 1.5, -2.3, 1.25;
     -2, 1.5, 30, 2, 2;
     4, -3.3, 2, 20, -1;
     -3.25, 1.25, 2, -1, 13];

b = [1; -1; 2; -3; 2.5];

n = length(b);
x = zeros(n, 1);
epsilon = 1e-7;
max_iter = 1000;

% Метод Зейделя
for k = 1:max_iter
    x_old = x;
    for i = 1:n
        sum1 = A(i, 1:i-1) * x(1:i-1);
        sum2 = A(i, i+1:n) * x_old(i+1:n);
        x(i) = (b(i) - sum1 - sum2) / A(i, i);
    end
    
    if norm(x - x_old, inf) < epsilon
        break;
    end
end

disp(x)

% Перевірка методом LU-розкладу
[L, U, P] = lu(A);
y = L \ (P * b);
x_lu = U \ y;
disp("Перевірка")
disp(x_lu)


%  --- Задача #2 ---


function[x, i, err] = m_njut(f, df, x0, max_i, tol)
    x = x0;
    for i = 1:max_i
        x_new = x - f(x) / df(x);
        err = abs(x_new - x); % різниця між новим і попереднім значенням
        x = x_new; % оновлюємо значення
        if err < tol
            return;
        end
    end
end

f = @(x) x^4 - 4*x^2 + 5*x - 20.1;
df = @(x) 4*x^3 - 8*x + 5;

x0 = 2;
max_i = 100;
tol = 1e-7;

% Виклик методу Ньютона
[x_njut, iter, err] = m_njut(f, df, x0, max_i, tol);

disp(["Метод Ньютона: x = ", num2str(x_njut)])
disp(["Ітерацій: ", num2str(iter)])
disp(["Похибка: ", num2str(err)])

% Перевірка за допомогою функції fzero
x_fzero = fzero(f, x0);
disp(["Результат fzero: x = ", num2str(x_fzero)]);

% Перевірка за допомогою функції fminbnd (пошук мінімуму в межах)
x_fminbnd = fminbnd(f, -10, 10);
disp(["Результат fminbnd: x = ", num2str(x_fminbnd)]);


%  --- Задача #3 ---


a = 0.4;
b = 0.4;
h = 0.1;
x = 2:h:3;

f = @(x) a + b ./ x;

% Інтерполяція функції
p = polyfit(x, f(x), length(x)-1); 
disp("Коефіцієнти інтерполяційного полінома: ");
disp(p);

% Права різницева похідна для середньої точки
mid_index = floor(length(x) / 2);
x_mid = x(mid_index); 
f_prime_right = (f(x(mid_index + 1)) - f(x_mid)) / h; 
disp(["Права різницева похідна в точці x = ", num2str(x_mid), ": ", num2str(f_prime_right)]);

% Наближене обчислення інтегралу методом правих прямокутників
integral_value = 0;
for i = 2:length(x)
    integral_value += f(x(i)) * h;
end
disp(["Наближене значення інтегралу: ", num2str(integral_value)]);


%  --- Задача #4 ---


a = 0.4;
b = 0.4;
h = 0.1;
t_start = 1;
t_end = 5;
x1 = 1.1;

f = @(x) a + b ./ x;

% Визначення правої частини диференціального рівняння
f_prime = @(t, x) 2.35 * sin(t * f(x));

% Явний метод Ейлера
t = t_start:h:t_end;
x = zeros(1, length(t));
x(1) = x1;

for i = 1:(length(t)-1)
    x(i+1) = x(i) + h * f_prime(t(i), x(i));
end

disp("Значення x(t) для кожного t:");
disp(x);

% Побудова графіка
plot(t, x, "-o");
xlabel("t");
ylabel("x(t)");
title("Наближений розв’язок задачі Коші методом Ейлера");
grid on;
