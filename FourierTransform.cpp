#include <vector>
#include <string>
#include "math.h"
#include <iostream>
#include "AudioFile.h"
#include <chrono>

using namespace std::chrono;

using std::vector, std::string;

// my own implementation of complex numbers
class Complex
{
public:
    double re;
    double im;

    Complex()
    {
        re = 0;
        im = 0;
    }

    Complex(double re, double im)
    {
        this->re = re;
        this->im = im;
    }

    Complex operator+(Complex const &obj)
    {
        return Complex(re + obj.re, im + obj.im);
    }

    Complex operator-(Complex const &obj)
    {
        return Complex(re - obj.re, im - obj.im);
    }

    Complex operator*(Complex const &obj)
    {
        return Complex(re * obj.re - im * obj.im, im * obj.re + re * obj.im);
    }

    void operator+=(Complex const &obj)
    {
        re += obj.re;
        im += obj.im;
    }
};

const double PI = 3.14159265358979;
void PrintVector(const vector<double> &vec)
{
    for (const double &x : vec)
    {
        std::cout << x << " ";
    }

    std::cout << "\n";
}

void PrintVector(const vector<Complex> &vec)
{
    for (const Complex &c : vec)
    {
        std::cout << "(" << c.re << "," << c.im << ") ";
    }

    std::cout << "\n";
}

class SlowFourierTransform
{
public:
    vector<Complex> coefficients;

    SlowFourierTransform(const vector<double> &values)
    {
        size_t N = values.size();
        coefficients.reserve(N);

        for (size_t k = 0; k < N; k++)
        {
            Complex x;

            for (size_t n = 0; n < N; n++)
            {
                double angle = 2 * k * n * PI / N;
                x += Complex(values[n] * cos(angle), values[n] * -sin(angle));
            }

            coefficients.push_back(x);
        }
    }

    vector<double> InverseDFT()
    {
        size_t N = coefficients.size();
        vector<double> values(N, 0);

        for (size_t k = 0; k < N; k++)
        {
            double x = 0;

            for (size_t n = 0; n < N; n++)
            {
                double angle = 2 * k * n * PI / N;
                x += (coefficients[n] * Complex(cos(angle), sin(angle))).re;
            }

            values[k] = x / N;
        }

        return values;
    }
};

std::pair<vector<Complex>, vector<Complex>> even_odd(vector<Complex> values)
{
    size_t N = values.size();
    vector<Complex> even;
    vector<Complex> odd;
    even.reserve(N / 2);
    odd.reserve(N / 2);

    for (size_t i = 0; i < N; i++)
    {
        if (i % 2 == 0)
        {
            even.push_back(values[i]);
        }
        else
        {
            odd.push_back(values[i]);
        }
    }

    return std::pair(even, odd);
}
class FastFourierTransform
{
public:
    vector<Complex> coefs;
    size_t initial_size;

    FastFourierTransform(const vector<double> &init_values)
    {
        // zero padding to be power of 2
        vector<double> values = init_values;
        size_t initial_n = init_values.size();
        initial_size = initial_n;

        size_t N = 1;
        while (N < initial_n)
        {
            N <<= 1;
        }
        coefs = vector<Complex>(N, Complex());
        for (size_t i = 0; i < initial_n; i++)
        {
            coefs[i].re = init_values[i];
        }
        coefs = FFT(coefs, N);
    }

    vector<Complex> FFT(vector<Complex> values, size_t n, bool isInverse = false)
    {
        if (n == 1)
        {
            return values;
        }
        else
        {
            auto [firstHalf, secondHalf] = even_odd(values);
            firstHalf = FFT(firstHalf, n / 2, isInverse);
            secondHalf = FFT(secondHalf, n / 2, isInverse);

            vector<Complex> ret;
            ret.resize(n);
            for (size_t k = 0; k < n / 2; k++)
            {
                Complex p = firstHalf[k];
                double angle = 2 * PI * k / n;
                Complex q = secondHalf[k] * Complex(cos(angle), (isInverse ? 1 : -1) * sin(angle));
                ret[k] = p + q;
                ret[k + (n / 2)] = p - q;
            }

            return ret;
        }
    }

    vector<double> InverseFFT()
    {
        size_t N = coefs.size();
        vector<Complex> inverseCoefs = FFT(coefs, N, true);
        vector<double> realCoefs;
        realCoefs.reserve(initial_size);
        for (size_t i = 0; i < initial_size; i++)
        {
            if (i < N)
            {
                realCoefs.push_back(inverseCoefs[i].re / N);
            }
            else
            {
                realCoefs.push_back(0);
            }
        }

        return realCoefs;
    }
};

// remove all frequencies from the audio file that are less that minFrequency or greater than maxFrequency
// This is an implementation of a bandpass filter
void Filter(AudioFile<double> &audio, double minFrequency, double maxFrequency)
{
    double seconds = audio.getLengthInSeconds();
    FastFourierTransform audio_fft(audio.samples[0]);

    for (size_t i = 0; i < minFrequency * seconds; i++)
    {
        audio_fft.coefs[i] = Complex();
    }
    for (size_t i = maxFrequency * seconds; i < audio_fft.coefs.size(); i++)
    {
        audio_fft.coefs[i] = Complex();
    }
    vector<double> filteredCoefs = audio_fft.InverseFFT();
    vector<vector<double>> newBuffer;
    newBuffer.push_back(filteredCoefs);
    bool ok = audio.setAudioBuffer(newBuffer);

    if (!ok)
    {
        throw("Error setting audio file buffer");
    }
}
int main(int argc, char **args)
{

    AudioFile<double> audioFile;
    audioFile.load("__audio__file__name");
    Filter(audioFile, 50, 2000); // filter audio file to only include frequencies between 50hz and 2000hz
    audioFile.save("__new__audio__file__name");
}