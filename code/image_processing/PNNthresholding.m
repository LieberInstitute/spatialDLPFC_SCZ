 net = denoisingNetwork('DnCNN'); %Load the pretrained denoising convolutional neural network, 'DnCNN'.
 dWFA = denoiseImage(WFA,net); %denoise the image