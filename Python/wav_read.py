import matplotlib.pyplot as plt
import numpy as np
import wave

file = 'input.wav'

with wave.open(file,'r') as wav:

    signal = wav.readframes(-1)
    signal = np.frombuffer(signal, np.int16)

    channels = [[] for channel in range(wav.getnchannels())]
    for index, datum in enumerate(signal):
        channels[index%len(channels)].append(datum)

    frames = wav.getframerate()
    Time=np.linspace(0, len(signal)/len(channels)/frames, num=int(len(signal)/len(channels)))

    plt.figure(1)
    plt.title('Signal in wav file')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    for channel in channels:
        plt.plot(Time,channel)
    plt.show()
