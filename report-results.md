# Lifting Factorization of PyWavelets filter banks

**Andrii Kryvyi** @ *Ukrainian Catholic University*

**Correspondence:** kryvy.pn@ucu.edu.ua  <br/>
**Date**: March 2026

### Summary

We have successfuly factored PyWavelets `float64`-valued filter banks into lifting steps over $\mathbb{Q}$ field. We have used algorithm of Sweldens et al. with several modifications. This factorization was done by a fully automatized pipeline built on top of SageMath computer algebra system.

## 1. Notes

To keep this report simple we avoid describing concepts of laurent polynomials and 2x2 polyphase matrices.

## 2. Pipeline input and output

Our pipeline takes as an input a filter bank consisting of two deconstruction FIR filters with polyphase matrix $P$ s.t. $\det P \approx 1$. This condition is satisfied by all discrete wavelet filters from PyWavelets catalog.

On the output it generates two integer delay terms (even and odd) and a sequence of lifting steps of 5 types parametrized by polynomial $q$. In `JSON` file step is stored as follows:

**Q-field format**
```json
{
  "type": "predict",
  "shift": 0,
  "coefficients": [
    {
      "numerator": 1,
      "denominator": 1
    }
  ]
}
```

**FP64 format**
```json
{
  "type": "predict",
  "shift": 0,
  "coefficients": [
    1
  ]
}
```

So polynomial $q$ can be recovered as follows:

$$
q = \sum_{i=0}^{\text{len}(\text{coefficients}) - 1} z^{i + \text{shift}} \cdot \text{coefficients}[i]
$$

Those steps can be described as right-multiplication by a matrix $Q$ that depends on polynomial $q$:

- **Predict**: predict odd from even

$$
Q_i = \begin{bmatrix}1 & q_i \\
0 & 1\end{bmatrix}
$$

- **Update**: update even from odd

$$
Q_i = \begin{bmatrix}1 & 0 \\
q_i & 1\end{bmatrix}
$$

- **Scale Even**: scale even by a constant factor (note: $q$ is a zero-degree constant polynomial)

$$
Q_i = \begin{bmatrix}q & 0 \\
0 & 1\end{bmatrix}
$$

- **Scale Odd**: scale odd by a constant factor (note: $q$ is a zero-degree constant polynomial)

$$
Q_i = \begin{bmatrix}1 & 0 \\
0 & q\end{bmatrix}
$$

- **Swap**: swap even and odd (note: no dependency on q)

$$
Q_i = \begin{bmatrix}0 & 1 \\
1 & 0\end{bmatrix}
$$

Thus one-level wavelet transform can be represented by repeated vector-matrix multiplication with $s_e$ / $s_o$ being even/odd part transformed into Laurent polynomial and $a$ / $d$ approximation/details (low-frequency/high-frequency)

$$
\begin{bmatrix} a & d \end{bmatrix} = \begin{bmatrix} s_e & s_o \end{bmatrix} \prod_{i=1}^n Q_i
$$

And reconstruction given as:
$$
\begin{bmatrix} s_e & s_o \end{bmatrix} = \begin{bmatrix} a & d \end{bmatrix} \prod_{i=n}^1 Q_i^{-1}
$$

With $Q_i^{-1}$ being an inverse of $Q_i$ which is trivial to retrieve.

## 3. Factoring results

We have successfully factorized all 106/106 FIR filter banks present in PyWavelets catalog. We have measured L2 error between original polyphase matrix $P$ and resulting matrices $P'_\mathbb{Q}$ (when multiplying in $\mathbb{Q}$-field) and $P'_{\text{FP64}}$ (when multiplying using `float64`):

$$
P' = \prod_{i=1}^nQ_i
$$

And L2 error measured as follows from the residuals matrix $R = P' - P$:

$$
L_2 = \sqrt{\|R_{11}\|^2_2+\|R_{12}\|^2_2+\|R_{21}\|^2_2+\|R_{22}\|^2_2} \\
$$

With L2-norm defined for Laurent polynomials as L2-norm over vector of coefficients.

The following error-rate of factorization was achieved:

![Error histogram](./figs/l2_histogram.svg)

| Wavelet | $\mathbb{Q}$ | **FP64** | Wavelet | $\mathbb{Q}$ | **FP64** | Wavelet | $\mathbb{Q}$ | **FP64** |
|---|---|---|---|---|---|---|---|---|
| bior1.1 | <span style='color: green'>5.75e-17</span> | <span style='color: green'>1.92e-16</span> | bior1.3 | <span style='color: green'>5.80e-17</span> | <span style='color: green'>1.92e-16</span> | bior1.5 | <span style='color: green'>5.84e-17</span> | <span style='color: green'>1.94e-16</span> |
| bior2.2 | <span style='color: green'>3.50e-17</span> | <span style='color: green'>8.78e-17</span> | bior2.4 | <span style='color: green'>4.28e-17</span> | <span style='color: green'>1.42e-16</span> | bior2.6 | <span style='color: green'>1.20e-16</span> | <span style='color: green'>1.77e-16</span> |
| bior2.8 | <span style='color: green'>1.13e-16</span> | <span style='color: green'>1.80e-16</span> | bior3.1 | <span style='color: green'>9.97e-17</span> | <span style='color: green'>2.60e-16</span> | bior3.3 | <span style='color: green'>8.47e-17</span> | <span style='color: green'>2.98e-16</span> |
| bior3.5 | <span style='color: green'>2.20e-16</span> | <span style='color: green'>8.46e-16</span> | bior3.7 | <span style='color: green'>1.91e-16</span> | <span style='color: green'>3.10e-16</span> | bior3.9 | <span style='color: green'>2.48e-16</span> | <span style='color: green'>3.35e-16</span> |
| bior4.4 | <span style='color: yellow'>1.07e-11</span> | <span style='color: yellow'>1.07e-11</span> | bior5.5 | <span style='color: yellow'>1.01e-11</span> | <span style='color: yellow'>1.01e-11</span> | bior6.8 | <span style='color: yellow'>7.63e-13</span> | <span style='color: yellow'>7.63e-13</span> |
| coif1 | <span style='color: green'>1.56e-16</span> | <span style='color: green'>5.09e-16</span> | coif2 | <span style='color: green'>1.25e-16</span> | <span style='color: green'>3.64e-15</span> | coif3 | <span style='color: green'>2.51e-16</span> | <span style='color: yellow'>1.08e-14</span> |
| coif4 | <span style='color: green'>1.90e-16</span> | <span style='color: yellow'>3.88e-14</span> | coif5 | <span style='color: green'>1.10e-15</span> | <span style='color: yellow'>1.58e-13</span> | coif6 | <span style='color: green'>2.87e-15</span> | <span style='color: yellow'>8.12e-13</span> |
| coif7 | <span style='color: yellow'>2.55e-14</span> | <span style='color: yellow'>3.89e-12</span> | coif8 | <span style='color: yellow'>4.76e-14</span> | <span style='color: yellow'>1.65e-11</span> | coif9 | <span style='color: yellow'>1.59e-13</span> | <span style='color: yellow'>3.53e-11</span> |
| coif10 | <span style='color: yellow'>4.94e-13</span> | <span style='color: yellow'>7.40e-10</span> | coif11 | <span style='color: yellow'>1.71e-12</span> | <span style='color: yellow'>1.78e-09</span> | coif12 | <span style='color: yellow'>3.74e-12</span> | <span style='color: yellow'>5.09e-09</span> |
| coif13 | <span style='color: yellow'>3.65e-11</span> | <span style='color: yellow'>3.42e-08</span> | coif14 | <span style='color: yellow'>1.12e-10</span> | <span style='color: red'>1.99e-07</span> | coif15 | <span style='color: yellow'>2.41e-10</span> | <span style='color: red'>1.97e-07</span> |
| coif16 | <span style='color: yellow'>9.61e-10</span> | <span style='color: red'>1.48e-06</span> | coif17 | <span style='color: yellow'>1.83e-09</span> | <span style='color: red'>2.91e-06</span> | db1 | <span style='color: green'>5.75e-17</span> | <span style='color: green'>1.92e-16</span> |
| db2 | <span style='color: green'>1.37e-16</span> | <span style='color: green'>5.70e-16</span> | db3 | <span style='color: green'>1.30e-16</span> | <span style='color: green'>1.24e-15</span> | db4 | <span style='color: green'>1.21e-16</span> | <span style='color: green'>1.92e-15</span> |
| db5 | <span style='color: green'>7.51e-16</span> | <span style='color: yellow'>1.07e-14</span> | db6 | <span style='color: green'>1.26e-15</span> | <span style='color: green'>8.40e-15</span> | db7 | <span style='color: green'>6.98e-15</span> | <span style='color: yellow'>1.76e-14</span> |
| db8 | <span style='color: yellow'>7.19e-14</span> | <span style='color: yellow'>7.97e-14</span> | db9 | <span style='color: yellow'>3.60e-12</span> | <span style='color: yellow'>3.60e-12</span> | db10 | <span style='color: yellow'>1.37e-11</span> | <span style='color: yellow'>1.37e-11</span> |
| db11 | <span style='color: yellow'>3.37e-10</span> | <span style='color: yellow'>3.37e-10</span> | db12 | <span style='color: yellow'>1.73e-09</span> | <span style='color: yellow'>1.73e-09</span> | db13 | <span style='color: red'>1.12e-06</span> | <span style='color: red'>1.12e-06</span> |
| db14 | <span style='color: red'>1.07e-04</span> | <span style='color: red'>1.07e-04</span> | db15 | <span style='color: red'>5.58e-03</span> | <span style='color: red'>5.58e-03</span> | db16 | <span style='color: red'>3.75e-02</span> | <span style='color: red'>3.75e-02</span> |
| db17 | <span style='color: red'>1.76e-01</span> | <span style='color: red'>1.76e-01</span> | db18 | <span style='color: red'>1.68e-01</span> | <span style='color: red'>1.68e-01</span> | db19 | <span style='color: red'>1.59e-01</span> | <span style='color: red'>1.59e-01</span> |
| db20 | <span style='color: red'>1.51e-01</span> | <span style='color: red'>1.51e-01</span> | db21 | <span style='color: red'>1.44e-01</span> | <span style='color: red'>1.44e-01</span> | db22 | <span style='color: red'>1.38e-01</span> | <span style='color: red'>1.38e-01</span> |
| db23 | <span style='color: red'>1.32e-01</span> | <span style='color: red'>1.32e-01</span> | db24 | <span style='color: red'>1.26e-01</span> | <span style='color: red'>1.26e-01</span> | db25 | <span style='color: red'>1.21e-01</span> | <span style='color: red'>1.21e-01</span> |
| db26 | <span style='color: red'>1.15e-01</span> | <span style='color: red'>1.15e-01</span> | db27 | <span style='color: red'>9.02e-02</span> | <span style='color: red'>9.02e-02</span> | db28 | <span style='color: red'>2.19e-01</span> | <span style='color: red'>2.19e-01</span> |
| db29 | <span style='color: red'>2.32e-01</span> | <span style='color: red'>2.32e-01</span> | db30 | <span style='color: red'>2.11e-01</span> | <span style='color: red'>2.11e-01</span> | db31 | <span style='color: red'>2.04e-01</span> | <span style='color: red'>2.04e-01</span> |
| db32 | <span style='color: red'>1.97e-01</span> | <span style='color: red'>1.97e-01</span> | db33 | <span style='color: red'>1.91e-01</span> | <span style='color: red'>1.91e-01</span> | db34 | <span style='color: red'>1.85e-01</span> | <span style='color: red'>1.85e-01</span> |
| db35 | <span style='color: red'>1.80e-01</span> | <span style='color: red'>1.80e-01</span> | db36 | <span style='color: red'>1.73e-01</span> | <span style='color: red'>1.73e-01</span> | db37 | <span style='color: red'>4.49e-01</span> | <span style='color: red'>4.49e-01</span> |
| db38 | <span style='color: red'>2.42e-01</span> | <span style='color: red'>2.42e-01</span> | dmey | <span style='color: red'>6.64e-01</span> | <span style='color: red'>6.64e-01</span> | haar | <span style='color: green'>5.75e-17</span> | <span style='color: green'>1.92e-16</span> |
| rbio1.1 | <span style='color: green'>5.75e-17</span> | <span style='color: green'>1.92e-16</span> | rbio1.3 | <span style='color: green'>5.80e-17</span> | <span style='color: green'>1.92e-16</span> | rbio1.5 | <span style='color: green'>5.84e-17</span> | <span style='color: green'>1.94e-16</span> |
| rbio2.2 | <span style='color: green'>3.11e-17</span> | <span style='color: green'>2.22e-16</span> | rbio2.4 | <span style='color: green'>3.90e-17</span> | <span style='color: green'>1.36e-16</span> | rbio2.6 | <span style='color: green'>9.92e-17</span> | <span style='color: green'>1.42e-16</span> |
| rbio2.8 | <span style='color: green'>1.49e-16</span> | <span style='color: green'>1.62e-16</span> | rbio3.1 | <span style='color: green'>9.97e-17</span> | <span style='color: green'>2.73e-16</span> | rbio3.3 | <span style='color: green'>8.47e-17</span> | <span style='color: green'>2.98e-16</span> |
| rbio3.5 | <span style='color: green'>2.20e-16</span> | <span style='color: green'>8.46e-16</span> | rbio3.7 | <span style='color: green'>1.91e-16</span> | <span style='color: green'>3.10e-16</span> | rbio3.9 | <span style='color: green'>2.48e-16</span> | <span style='color: green'>3.35e-16</span> |
| rbio4.4 | <span style='color: yellow'>2.74e-12</span> | <span style='color: yellow'>2.74e-12</span> | rbio5.5 | <span style='color: yellow'>1.09e-08</span> | <span style='color: yellow'>1.09e-08</span> | rbio6.8 | <span style='color: yellow'>3.81e-13</span> | <span style='color: yellow'>3.81e-13</span> |
| sym2 | <span style='color: yellow'>8.77e-13</span> | <span style='color: yellow'>8.77e-13</span> | sym3 | <span style='color: yellow'>1.60e-11</span> | <span style='color: yellow'>1.60e-11</span> | sym4 | <span style='color: yellow'>1.95e-12</span> | <span style='color: yellow'>1.95e-12</span> |
| sym5 | <span style='color: yellow'>7.35e-13</span> | <span style='color: yellow'>7.35e-13</span> | sym6 | <span style='color: yellow'>5.32e-12</span> | <span style='color: yellow'>5.32e-12</span> | sym7 | <span style='color: yellow'>6.93e-10</span> | <span style='color: yellow'>6.93e-10</span> |
| sym8 | <span style='color: yellow'>3.10e-11</span> | <span style='color: yellow'>3.10e-11</span> | sym9 | <span style='color: yellow'>1.87e-11</span> | <span style='color: yellow'>1.87e-11</span> | sym10 | <span style='color: yellow'>5.66e-12</span> | <span style='color: yellow'>5.66e-12</span> |
| sym11 | <span style='color: yellow'>2.59e-10</span> | <span style='color: yellow'>2.59e-10</span> | sym12 | <span style='color: yellow'>2.20e-12</span> | <span style='color: yellow'>2.20e-12</span> | sym13 | <span style='color: yellow'>8.37e-09</span> | <span style='color: yellow'>8.37e-09</span> |
| sym14 | <span style='color: yellow'>1.78e-12</span> | <span style='color: yellow'>1.78e-12</span> | sym15 | <span style='color: yellow'>1.65e-11</span> | <span style='color: yellow'>1.65e-11</span> | sym16 | <span style='color: yellow'>2.87e-11</span> | <span style='color: yellow'>2.87e-11</span> |
| sym17 | <span style='color: yellow'>4.58e-09</span> | <span style='color: yellow'>4.58e-09</span> | sym18 | <span style='color: yellow'>2.18e-10</span> | <span style='color: yellow'>2.18e-10</span> | sym19 | <span style='color: yellow'>4.16e-11</span> | <span style='color: yellow'>4.16e-11</span> |
| sym20 | <span style='color: yellow'>3.05e-10</span> | <span style='color: yellow'>3.05e-10</span> |

The above table is colored as: `green` if $L_2 < 10^{-14}$, `yellow` if $L_2 < 10^{-7}$ and otherwise `red`.

![Error scatter](./figs/l2_scatter.svg)

We can see that FP64 error is highly correlated with Q-field error, thus we can conclude that factorization is stable in terms of recovering original polyphase matrix, because almost all error comes from projecting input matrix into space of matrices with determinant exactly *1*, but not from converting all coefficient to `float64`.

![Error scatter](./figs/l2_scatter_coif.svg)

Only for coiflets the datatype matters and converting fractions to `float64` increases error rate by a constant factor.

## 4. Comparing numerical stability


