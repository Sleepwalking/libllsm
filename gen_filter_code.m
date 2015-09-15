#!/usr/bin/octave

step_freq = 0.05;
number = 18;
order = 5;
ripple = 1.0;
coef_size = order + 1;

filter_l_a = [];
filter_l_b = [];
filter_h_a = [];
filter_h_b = [];

for i = 1:number
  fc = i * step_freq;
  [b a] = cheby1(order, ripple, fc);
  filter_l_a = [filter_l_a a];
  filter_l_b = [filter_l_b b];
  [b a] = cheby1(order, ripple, fc, 'high');
  filter_h_a = [filter_h_a a];
  filter_h_b = [filter_h_b b];
end

function printfarr(arr, column)
  n = length(arr) / column;
  for i = 1:n
    printf('  %11.9f, ', arr((i - 1) * column + 1));
    for j = 2:column - 1
      printf('%13.9f, ', arr((i - 1) * column + j));
    end
    if(i == n)
      printf('%13.9f\n', arr(i * column));
    else
      printf('%13.9f,\n', arr(i * column));
    end
  end
end

printf('// The following are program-generated codes and are not supposed to be modified by hand.\n\n');
printf('#ifndef LLSM_FILTER_COEF\n');
printf('#define LLSM_FILTER_COEF\n\n');
printf('#include \"common.h\"\n\n');

printf('static const int coef_size = %d;\n', coef_size);
printf('static const int filter_number = %d;\n', number);
printf('static const FP_TYPE step_freq = %f;\n\n', step_freq);

printf('static const FP_TYPE cheby_l_a[%d] = {\n', length(filter_l_a));
printfarr(filter_l_a, coef_size);
printf('};\n\n');

printf('static const FP_TYPE cheby_l_b[%d] = {\n', length(filter_l_b));
printfarr(filter_l_b, coef_size);
printf('};\n\n');

printf('static const FP_TYPE cheby_h_a[%d] = {\n', length(filter_h_a));
printfarr(filter_h_a, coef_size);
printf('};\n\n');

printf('static const FP_TYPE cheby_h_b[%d] = {\n', length(filter_h_b));
printfarr(filter_h_b, coef_size);
printf('};\n\n');

printf('#endif\n\n');

