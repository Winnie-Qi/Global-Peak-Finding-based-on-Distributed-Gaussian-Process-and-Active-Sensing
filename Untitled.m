x_train = [-4, -4]; % ѵ��������
x_test = [-3.5, -3.5]; % ���Ե�����

% �����˹���̵Ĳ���
E = 10; % KL Kernel Expansion��ά��
sigma_v = 0.1; % ��˹���̵���������

% ����ѵ����ĺ˾��� G
G = zeros(size(x_train, 1), E);
for i = 1:E
    for j = 1:E
        phi = 2 * sin((i-0.5)*pi*x_train(:,1)) .* sin((j-0.5)*pi*x_train(:,2));
        lambda = 16 / (((2*i-1)^2*pi^2) * ((2*j-1)^2*pi^2));
        G(:, (i-1)*E+j) = phi * sqrt(lambda);
    end
end

% ������Ե�ĺ˾��� Phi
Phi = zeros(size(x_test, 1), E);
for i = 1:E
    for j = 1:E
        phi = 2 * sin((i-0.5)*pi*x_test(:,1)) .* sin((j-0.5)*pi*x_test(:,2));
        Phi(:, (i-1)*E+j) = phi;
    end
end

% ���� H_E
N = size(x_train, 1);
m = size(x_train, 2);
Lambda_E = diag(1 ./ (((2*(1:E)-1).^2*pi^2).^2));
H_E = (G' * G / (N*m) + sigma_v^2 / (N*m) * Lambda_E)^(-1) * G' / (N*m);

% ����ѵ����ĸ߶����� y
y_train = f(-3.5,-3.5);

% ������Ե�ĸ߶����� f_E(x')
f_E = Phi * H_E * y_train;