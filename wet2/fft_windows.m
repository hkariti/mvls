function [DFT_window] =  fft_windows( x , N)
windowed= im2col(x,[1 N],'distinct');
DFT_window=fft(windowed);
DFT_window= DFT_window';
end