syms y R
f = -1/2 * (tanh(R*(y-1/3))+tanh(R*(y+1/3))); % Eaves2016JFM
taylor_series = taylor(f, y, 'Order', 4); % Taylor series expansion up to 4th order
disp(taylor_series)