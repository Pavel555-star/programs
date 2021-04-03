import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wave

file = 'input.wav'

with wave.open(file,'r') as wav:
    signal = wav.readframes(-1)
    signal = np.frombuffer(signal, np.int16)

    channels = [[] for channel in range(wav.getnchannels())]
    for index, data in enumerate(signal):
        channels[index%len(channels)].append(data)

    frames = wav.getframerate()
    Time=np.linspace(0, len(signal)/len(channels)/frames, num=int(len(signal)/len(channels)))

    plt.figure(1)
    plt.title('Fourier transform')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    i = 0
    for channel in channels:
        y = np.fft.fft(channel)
        x = np.fft.fftfreq(y.shape[-1]) * frames
        plt.plot(x, y.imag, x, y.real)
        
        csv = {'Frequency (Hz)': x, 'Amplitude real': y.real, 'Amplitude imaginary': y.imag}
        output = pd.DataFrame(csv)
        output.to_csv('output' + str(i) + '.csv')
        i = i + 1
    plt.show()
