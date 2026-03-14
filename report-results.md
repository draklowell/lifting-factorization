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

We have successfully factorized all 106/106 FIR filter banks present in PyWavelets catalog. We have measured L2 error between original polyphase matrix $P$ and resulting matrix $P'$:

$$
P' = \prod_{i=1}^nQ_i
$$

And L2 error measured as follows from the residuals matrix $R = P' - P$:

$$
L_2 = \sqrt{\|R_{11}\|^2_2+\|R_{12}\|^2_2+\|R_{21}\|^2_2+\|R_{22}\|^2_2} \\
$$

With L2-norm defined for Laurent polynomials as L2-norm over vector of coefficients.

We consider two $P'$ reconstructed matrices: first by multiplication in $\mathbb{Q}$ field, and second by multiplication after converting all coefficients to `float64`. The following error-rate of factorization was achieved:

![Error histogram](./figs/l2_histogram.svg)

| Wavelet | $\mathbb{Q}$ | **FP64** | Wavelet | $\mathbb{Q}$ | **FP64** | Wavelet | $\mathbb{Q}$ | **FP64** |
|---|---|---|---|---|---|---|---|---|
| bior1.1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ | bior1.3 | $${\color{green} \text{5.80e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ | bior1.5 | $${\color{green} \text{5.84e-17}}$$ | $${\color{green} \text{1.94e-16}}$$ |
| bior2.2 | $${\color{green} \text{3.50e-17}}$$ | $${\color{green} \text{8.78e-17}}$$ | bior2.4 | $${\color{green} \text{4.28e-17}}$$ | $${\color{green} \text{1.42e-16}}$$ | bior2.6 | $${\color{green} \text{1.20e-16}}$$ | $${\color{green} \text{1.77e-16}}$$ |
| bior2.8 | $${\color{green} \text{1.13e-16}}$$ | $${\color{green} \text{1.80e-16}}$$ | bior3.1 | $${\color{green} \text{9.97e-17}}$$ | $${\color{green} \text{2.60e-16}}$$ | bior3.3 | $${\color{green} \text{8.47e-17}}$$ | $${\color{green} \text{2.98e-16}}$$ |
| bior3.5 | $${\color{green} \text{2.20e-16}}$$ | $${\color{green} \text{8.46e-16}}$$ | bior3.7 | $${\color{green} \text{1.91e-16}}$$ | $${\color{green} \text{3.10e-16}}$$ | bior3.9 | $${\color{green} \text{2.48e-16}}$$ | $${\color{green} \text{3.35e-16}}$$ |
| bior4.4 | $${\color{yellow} \text{1.07e-11}}$$ | $${\color{yellow} \text{1.07e-11}}$$ | bior5.5 | $${\color{yellow} \text{1.01e-11}}$$ | $${\color{yellow} \text{1.01e-11}}$$ | bior6.8 | $${\color{yellow} \text{7.63e-13}}$$ | $${\color{yellow} \text{7.63e-13}}$$ |
| coif1 | $${\color{green} \text{1.56e-16}}$$ | $${\color{green} \text{5.09e-16}}$$ | coif2 | $${\color{green} \text{1.25e-16}}$$ | $${\color{green} \text{3.64e-15}}$$ | coif3 | $${\color{green} \text{2.51e-16}}$$ | $${\color{yellow} \text{1.08e-14}}$$ |
| coif4 | $${\color{green} \text{1.90e-16}}$$ | $${\color{yellow} \text{3.88e-14}}$$ | coif5 | $${\color{green} \text{1.10e-15}}$$ | $${\color{yellow} \text{1.58e-13}}$$ | coif6 | $${\color{green} \text{2.87e-15}}$$ | $${\color{yellow} \text{8.12e-13}}$$ |
| coif7 | $${\color{yellow} \text{2.55e-14}}$$ | $${\color{yellow} \text{3.89e-12}}$$ | coif8 | $${\color{yellow} \text{4.76e-14}}$$ | $${\color{yellow} \text{1.65e-11}}$$ | coif9 | $${\color{yellow} \text{1.59e-13}}$$ | $${\color{yellow} \text{3.53e-11}}$$ |
| coif10 | $${\color{yellow} \text{4.94e-13}}$$ | $${\color{yellow} \text{7.40e-10}}$$ | coif11 | $${\color{yellow} \text{1.71e-12}}$$ | $${\color{yellow} \text{1.78e-09}}$$ | coif12 | $${\color{yellow} \text{3.74e-12}}$$ | $${\color{yellow} \text{5.09e-09}}$$ |
| coif13 | $${\color{yellow} \text{3.65e-11}}$$ | $${\color{yellow} \text{3.42e-08}}$$ | coif14 | $${\color{yellow} \text{1.12e-10}}$$ | $${\color{red} \text{1.99e-07}}$$ | coif15 | $${\color{yellow} \text{2.41e-10}}$$ | $${\color{red} \text{1.97e-07}}$$ |
| coif16 | $${\color{yellow} \text{9.61e-10}}$$ | $${\color{red} \text{1.48e-06}}$$ | coif17 | $${\color{yellow} \text{1.83e-09}}$$ | $${\color{red} \text{2.91e-06}}$$ | db1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ |
| db2 | $${\color{green} \text{1.37e-16}}$$ | $${\color{green} \text{5.70e-16}}$$ | db3 | $${\color{green} \text{1.30e-16}}$$ | $${\color{green} \text{1.24e-15}}$$ | db4 | $${\color{green} \text{1.21e-16}}$$ | $${\color{green} \text{1.92e-15}}$$ |
| db5 | $${\color{green} \text{7.51e-16}}$$ | $${\color{yellow} \text{1.07e-14}}$$ | db6 | $${\color{green} \text{1.26e-15}}$$ | $${\color{green} \text{8.40e-15}}$$ | db7 | $${\color{green} \text{6.98e-15}}$$ | $${\color{yellow} \text{1.76e-14}}$$ |
| db8 | $${\color{yellow} \text{7.19e-14}}$$ | $${\color{yellow} \text{7.97e-14}}$$ | db9 | $${\color{yellow} \text{3.60e-12}}$$ | $${\color{yellow} \text{3.60e-12}}$$ | db10 | $${\color{yellow} \text{1.37e-11}}$$ | $${\color{yellow} \text{1.37e-11}}$$ |
| db11 | $${\color{yellow} \text{3.37e-10}}$$ | $${\color{yellow} \text{3.37e-10}}$$ | db12 | $${\color{yellow} \text{1.73e-09}}$$ | $${\color{yellow} \text{1.73e-09}}$$ | db13 | $${\color{red} \text{1.12e-06}}$$ | $${\color{red} \text{1.12e-06}}$$ |
| db14 | $${\color{red} \text{1.07e-04}}$$ | $${\color{red} \text{1.07e-04}}$$ | db15 | $${\color{red} \text{5.58e-03}}$$ | $${\color{red} \text{5.58e-03}}$$ | db16 | $${\color{red} \text{3.75e-02}}$$ | $${\color{red} \text{3.75e-02}}$$ |
| db17 | $${\color{red} \text{1.76e-01}}$$ | $${\color{red} \text{1.76e-01}}$$ | db18 | $${\color{red} \text{1.68e-01}}$$ | $${\color{red} \text{1.68e-01}}$$ | db19 | $${\color{red} \text{1.59e-01}}$$ | $${\color{red} \text{1.59e-01}}$$ |
| db20 | $${\color{red} \text{1.51e-01}}$$ | $${\color{red} \text{1.51e-01}}$$ | db21 | $${\color{red} \text{1.44e-01}}$$ | $${\color{red} \text{1.44e-01}}$$ | db22 | $${\color{red} \text{1.38e-01}}$$ | $${\color{red} \text{1.38e-01}}$$ |
| db23 | $${\color{red} \text{1.32e-01}}$$ | $${\color{red} \text{1.32e-01}}$$ | db24 | $${\color{red} \text{1.26e-01}}$$ | $${\color{red} \text{1.26e-01}}$$ | db25 | $${\color{red} \text{1.21e-01}}$$ | $${\color{red} \text{1.21e-01}}$$ |
| db26 | $${\color{red} \text{1.15e-01}}$$ | $${\color{red} \text{1.15e-01}}$$ | db27 | $${\color{red} \text{9.02e-02}}$$ | $${\color{red} \text{9.02e-02}}$$ | db28 | $${\color{red} \text{2.19e-01}}$$ | $${\color{red} \text{2.19e-01}}$$ |
| db29 | $${\color{red} \text{2.32e-01}}$$ | $${\color{red} \text{2.32e-01}}$$ | db30 | $${\color{red} \text{2.11e-01}}$$ | $${\color{red} \text{2.11e-01}}$$ | db31 | $${\color{red} \text{2.04e-01}}$$ | $${\color{red} \text{2.04e-01}}$$ |
| db32 | $${\color{red} \text{1.97e-01}}$$ | $${\color{red} \text{1.97e-01}}$$ | db33 | $${\color{red} \text{1.91e-01}}$$ | $${\color{red} \text{1.91e-01}}$$ | db34 | $${\color{red} \text{1.85e-01}}$$ | $${\color{red} \text{1.85e-01}}$$ |
| db35 | $${\color{red} \text{1.80e-01}}$$ | $${\color{red} \text{1.80e-01}}$$ | db36 | $${\color{red} \text{1.73e-01}}$$ | $${\color{red} \text{1.73e-01}}$$ | db37 | $${\color{red} \text{4.49e-01}}$$ | $${\color{red} \text{4.49e-01}}$$ |
| db38 | $${\color{red} \text{2.42e-01}}$$ | $${\color{red} \text{2.42e-01}}$$ | dmey | $${\color{red} \text{6.64e-01}}$$ | $${\color{red} \text{6.64e-01}}$$ | haar | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ |
| rbio1.1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ | rbio1.3 | $${\color{green} \text{5.80e-17}}$$ | $${\color{green} \text{1.92e-16}}$$ | rbio1.5 | $${\color{green} \text{5.84e-17}}$$ | $${\color{green} \text{1.94e-16}}$$ |
| rbio2.2 | $${\color{green} \text{3.11e-17}}$$ | $${\color{green} \text{2.22e-16}}$$ | rbio2.4 | $${\color{green} \text{3.90e-17}}$$ | $${\color{green} \text{1.36e-16}}$$ | rbio2.6 | $${\color{green} \text{9.92e-17}}$$ | $${\color{green} \text{1.42e-16}}$$ |
| rbio2.8 | $${\color{green} \text{1.49e-16}}$$ | $${\color{green} \text{1.62e-16}}$$ | rbio3.1 | $${\color{green} \text{9.97e-17}}$$ | $${\color{green} \text{2.73e-16}}$$ | rbio3.3 | $${\color{green} \text{8.47e-17}}$$ | $${\color{green} \text{2.98e-16}}$$ |
| rbio3.5 | $${\color{green} \text{2.20e-16}}$$ | $${\color{green} \text{8.46e-16}}$$ | rbio3.7 | $${\color{green} \text{1.91e-16}}$$ | $${\color{green} \text{3.10e-16}}$$ | rbio3.9 | $${\color{green} \text{2.48e-16}}$$ | $${\color{green} \text{3.35e-16}}$$ |
| rbio4.4 | $${\color{yellow} \text{2.74e-12}}$$ | $${\color{yellow} \text{2.74e-12}}$$ | rbio5.5 | $${\color{yellow} \text{1.09e-08}}$$ | $${\color{yellow} \text{1.09e-08}}$$ | rbio6.8 | $${\color{yellow} \text{3.81e-13}}$$ | $${\color{yellow} \text{3.81e-13}}$$ |
| sym2 | $${\color{yellow} \text{8.77e-13}}$$ | $${\color{yellow} \text{8.77e-13}}$$ | sym3 | $${\color{yellow} \text{1.60e-11}}$$ | $${\color{yellow} \text{1.60e-11}}$$ | sym4 | $${\color{yellow} \text{1.95e-12}}$$ | $${\color{yellow} \text{1.95e-12}}$$ |
| sym5 | $${\color{yellow} \text{7.35e-13}}$$ | $${\color{yellow} \text{7.35e-13}}$$ | sym6 | $${\color{yellow} \text{5.32e-12}}$$ | $${\color{yellow} \text{5.32e-12}}$$ | sym7 | $${\color{yellow} \text{6.93e-10}}$$ | $${\color{yellow} \text{6.93e-10}}$$ |
| sym8 | $${\color{yellow} \text{3.10e-11}}$$ | $${\color{yellow} \text{3.10e-11}}$$ | sym9 | $${\color{yellow} \text{1.87e-11}}$$ | $${\color{yellow} \text{1.87e-11}}$$ | sym10 | $${\color{yellow} \text{5.66e-12}}$$ | $${\color{yellow} \text{5.66e-12}}$$ |
| sym11 | $${\color{yellow} \text{2.59e-10}}$$ | $${\color{yellow} \text{2.59e-10}}$$ | sym12 | $${\color{yellow} \text{2.20e-12}}$$ | $${\color{yellow} \text{2.20e-12}}$$ | sym13 | $${\color{yellow} \text{8.37e-09}}$$ | $${\color{yellow} \text{8.37e-09}}$$ |
| sym14 | $${\color{yellow} \text{1.78e-12}}$$ | $${\color{yellow} \text{1.78e-12}}$$ | sym15 | $${\color{yellow} \text{1.65e-11}}$$ | $${\color{yellow} \text{1.65e-11}}$$ | sym16 | $${\color{yellow} \text{2.87e-11}}$$ | $${\color{yellow} \text{2.87e-11}}$$ |
| sym17 | $${\color{yellow} \text{4.58e-09}}$$ | $${\color{yellow} \text{4.58e-09}}$$ | sym18 | $${\color{yellow} \text{2.18e-10}}$$ | $${\color{yellow} \text{2.18e-10}}$$ | sym19 | $${\color{yellow} \text{4.16e-11}}$$ | $${\color{yellow} \text{4.16e-11}}$$ |
| sym20 | $${\color{yellow} \text{3.05e-10}}$$ | $${\color{yellow} \text{3.05e-10}}$$ |

The above table is colored as: `green` if $L_2 < 10^{-14}$, `yellow` if $L_2 < 10^{-7}$ and otherwise `red`.

![Error scatter](./figs/l2_scatter.svg)

We can see that FP64 error is highly correlated with Q-field error, thus we can conclude that factorization is stable in terms of recovering original polyphase matrix, because almost all error comes from projecting input matrix into space of matrices with determinant exactly *1*, but not from converting all coefficient to `float64`.

![Error scatter](./figs/l2_scatter_coif.svg)

Only for coiflets the datatype matters and converting fractions to `float64` increases error rate by a constant factor.

## 4. Comparing numerical stability


