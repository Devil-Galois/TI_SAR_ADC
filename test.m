clear;
clc;
close all;

try
    adc = adc_model();
    fprintf('ADC Initial Successfully!\n');
catch ME
    error('ADC Initial Failed\nMessage: %s', ME.message);
end

adc.model_init();
adc.run_analysis();
fprintf("========All End========");