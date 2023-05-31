clc; clear all; 
%начальные значения (параметр регуляризации альфа,число Куранта бета и кол-во ячеек N)
alphaGidro = 0.2; betaGidro = 0.2; NxGidro =400;
alphaGas = 0.1; betaGas = 0.2; NxGas = 400;

%гидроДин
gidroResult = gidro(alphaGidro, betaGidro, NxGidro);
xxGidro = gidroResult(1,:);
roGidro = gidroResult(2,:);

%средняя схема гидро дин
%gidroResultMiddle = gidro_middle(alphaGidro, betaGidro, NxGidro);
%xxGidroMiddle = gidroResultMiddle(1,:);
%roGidroMiddle = gidroResultMiddle(2,:);

%газ
gasResult = gas(alphaGas, betaGas, NxGas);  
xxGas = gasResult(1,:);
roGas = gasResult(2,:);

%среднее газ 
%gasResultMiddle = gas_middle(alphaGas, betaGas, NxGas);
%xxGasMiddle = gasResultMiddle(1,:);
%roGasMiddle = gasResultMiddle(2,:);

%уравнения состояния
%gidroNewEqResult = gidro_new_eq(alphaGidro, betaGidro, NxGidro);
%xxGidroNew = gidroNewResult(1,:);
%roGidroNew = gidroNewResult(2,:);

 
%графики
hold on;
plot(xxGidro,roGidro,'r', 'LineWidth',2); 
%plot(xxGidroMiddle,roGidroMiddle,'--b', 'LineWidth',2); 
plot(xxGas, roGas, '--b','LineWidth',2); 
%plot(xxGasMiddle, roGasMiddle, '--k');


%plot(xxGidro,roGidro,'r', 'LineWidth',2); 

xlabel('x')
ylabel('\rho')
%legend('гидро (направленная)','газ направленная)')
grid minor
hold off;
