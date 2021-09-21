function waveform = resample(path, fs, fs_res)
%RESAMPLE Summary of this function goes here
%   Detailed explanation goes here
    t_int = 0.1;    % [s]
    file = fopen(path);
    data = fread(file, 2*t_int*fs, 'int16');
    waveform_rx = data(1:2:end).' + 1j*data(2:2:end).';
    [up, down] = rat((1/fs_res) / (1/fs));
    waveform = resample(waveform_rx, down, up);
end

