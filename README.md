# Lifting Factorization of PyWavelets filter banks

**Andrii Kryvyi** @ *Ukrainian Catholic University*

**Correspondence:** kryvy.pn@ucu.edu.ua  <br/>
**Date**: March 2026

### Summary

We have successfully factored PyWavelets `float64`-valued filter banks into lifting steps over $\mathbb{Q}$ field. We have used algorithm of Sweldens et al. with several modifications. This factorization was done by a fully automatized pipeline built on top of SageMath computer algebra system.

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

We produced lifting decompositions for all 106/106 PyWavelets FIR filter banks, but the reconstruction accuracy varies substantially across families. The factorization is near-exact for many biorthogonal, reverse-biorthogonal, low-order Daubechies, and some coiflet/symlet filters, while higher-order Daubechies and dmey show large residuals.

To validate factorization, we have measured L2 error between original polyphase matrix $P$ and resulting matrix $P'$:

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
| bior1.1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ | bior1.3 | $${\color{green} \text{5.80e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ | bior1.5 | $${\color{green} \text{5.84e-17}}$$ | $${\color{green} \text{7.48e-17}}$$ |
| bior2.2 | $${\color{green} \text{3.50e-17}}$$ | $${\color{green} \text{4.98e-17}}$$ | bior2.4 | $${\color{green} \text{4.28e-17}}$$ | $${\color{green} \text{1.29e-16}}$$ | bior2.6 | $${\color{green} \text{1.20e-16}}$$ | $${\color{green} \text{2.15e-16}}$$ |
| bior2.8 | $${\color{green} \text{1.13e-16}}$$ | $${\color{green} \text{1.72e-16}}$$ | bior3.1 | $${\color{green} \text{9.97e-17}}$$ | $${\color{green} \text{2.26e-16}}$$ | bior3.3 | $${\color{green} \text{8.47e-17}}$$ | $${\color{green} \text{1.75e-16}}$$ |
| bior3.5 | $${\color{green} \text{2.20e-16}}$$ | $${\color{green} \text{1.31e-16}}$$ | bior3.7 | $${\color{green} \text{1.91e-16}}$$ | $${\color{green} \text{1.44e-16}}$$ | bior3.9 | $${\color{green} \text{2.48e-16}}$$ | $${\color{green} \text{2.58e-16}}$$ |
| bior4.4 | $${\color{yellow} \text{1.07e-11}}$$ | $${\color{yellow} \text{1.07e-11}}$$ | bior5.5 | $${\color{yellow} \text{1.01e-11}}$$ | $${\color{yellow} \text{1.01e-11}}$$ | bior6.8 | $${\color{yellow} \text{7.63e-13}}$$ | $${\color{yellow} \text{7.63e-13}}$$ |
| coif1 | $${\color{green} \text{1.15e-16}}$$ | $${\color{green} \text{1.85e-16}}$$ | coif10 | $${\color{green} \text{8.46e-17}}$$ | $${\color{green} \text{4.88e-16}}$$ | coif11 | $${\color{green} \text{1.60e-16}}$$ | $${\color{green} \text{6.28e-16}}$$ |
| coif12 | $${\color{green} \text{1.50e-16}}$$ | $${\color{green} \text{5.90e-16}}$$ | coif13 | $${\color{green} \text{1.19e-16}}$$ | $${\color{green} \text{1.91e-15}}$$ | coif14 | $${\color{green} \text{1.57e-16}}$$ | $${\color{green} \text{8.63e-16}}$$ |
| coif15 | $${\color{green} \text{1.38e-16}}$$ | $${\color{green} \text{1.68e-15}}$$ | coif16 | $${\color{green} \text{2.74e-16}}$$ | $${\color{green} \text{8.38e-16}}$$ | coif17 | $${\color{green} \text{7.37e-17}}$$ | $${\color{green} \text{3.16e-15}}$$ |
| coif2 | $${\color{green} \text{8.59e-17}}$$ | $${\color{green} \text{2.65e-16}}$$ | coif3 | $${\color{green} \text{1.19e-16}}$$ | $${\color{green} \text{2.32e-16}}$$ | coif4 | $${\color{green} \text{9.53e-17}}$$ | $${\color{green} \text{1.80e-16}}$$ |
| coif5 | $${\color{green} \text{6.23e-17}}$$ | $${\color{green} \text{2.97e-16}}$$ | coif6 | $${\color{green} \text{1.65e-16}}$$ | $${\color{green} \text{4.26e-16}}$$ | coif7 | $${\color{green} \text{1.53e-16}}$$ | $${\color{green} \text{4.16e-16}}$$ |
| coif8 | $${\color{green} \text{1.37e-16}}$$ | $${\color{green} \text{3.01e-16}}$$ | coif9 | $${\color{green} \text{6.40e-17}}$$ | $${\color{green} \text{8.77e-16}}$$ | db1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ |
| db10 | $${\color{green} \text{1.32e-16}}$$ | $${\color{green} \text{2.24e-16}}$$ | db11 | $${\color{green} \text{1.43e-16}}$$ | $${\color{green} \text{2.08e-16}}$$ | db12 | $${\color{green} \text{1.19e-16}}$$ | $${\color{green} \text{2.81e-16}}$$ |
| db13 | $${\color{green} \text{1.38e-16}}$$ | $${\color{green} \text{2.33e-16}}$$ | db14 | $${\color{green} \text{7.35e-17}}$$ | $${\color{green} \text{3.87e-16}}$$ | db15 | $${\color{green} \text{2.71e-17}}$$ | $${\color{green} \text{3.25e-16}}$$ |
| db16 | $${\color{green} \text{5.38e-17}}$$ | $${\color{green} \text{3.31e-16}}$$ | db17 | $${\color{green} \text{1.57e-16}}$$ | $${\color{green} \text{3.53e-16}}$$ | db18 | $${\color{green} \text{1.40e-16}}$$ | $${\color{green} \text{3.03e-16}}$$ |
| db19 | $${\color{green} \text{9.13e-17}}$$ | $${\color{green} \text{3.65e-16}}$$ | db2 | $${\color{green} \text{1.37e-16}}$$ | $${\color{green} \text{7.51e-17}}$$ | db20 | $${\color{green} \text{7.82e-17}}$$ | $${\color{green} \text{3.17e-16}}$$ |
| db21 | $${\color{green} \text{4.58e-17}}$$ | $${\color{green} \text{2.60e-16}}$$ | db22 | $${\color{green} \text{1.07e-16}}$$ | $${\color{green} \text{5.26e-16}}$$ | db23 | $${\color{green} \text{1.43e-16}}$$ | $${\color{green} \text{3.57e-16}}$$ |
| db24 | $${\color{green} \text{9.61e-17}}$$ | $${\color{green} \text{3.53e-16}}$$ | db25 | $${\color{green} \text{9.28e-17}}$$ | $${\color{green} \text{4.63e-16}}$$ | db26 | $${\color{green} \text{6.42e-17}}$$ | $${\color{green} \text{4.12e-16}}$$ |
| db27 | $${\color{green} \text{9.60e-17}}$$ | $${\color{green} \text{3.66e-16}}$$ | db28 | $${\color{green} \text{6.73e-17}}$$ | $${\color{green} \text{3.33e-16}}$$ | db29 | $${\color{green} \text{8.05e-17}}$$ | $${\color{green} \text{1.17e-15}}$$ |
| db3 | $${\color{green} \text{1.24e-16}}$$ | $${\color{green} \text{1.78e-16}}$$ | db30 | $${\color{green} \text{8.94e-17}}$$ | $${\color{green} \text{1.60e-15}}$$ | db31 | $${\color{green} \text{9.25e-17}}$$ | $${\color{green} \text{6.72e-16}}$$ |
| db32 | $${\color{green} \text{1.04e-16}}$$ | $${\color{green} \text{3.54e-16}}$$ | db33 | $${\color{green} \text{8.86e-17}}$$ | $${\color{green} \text{3.14e-16}}$$ | db34 | $${\color{green} \text{1.07e-16}}$$ | $${\color{green} \text{5.87e-16}}$$ |
| db35 | $${\color{green} \text{1.01e-16}}$$ | $${\color{green} \text{5.09e-16}}$$ | db36 | $${\color{green} \text{1.01e-16}}$$ | $${\color{green} \text{7.37e-16}}$$ | db37 | $${\color{green} \text{9.42e-17}}$$ | $${\color{green} \text{9.95e-16}}$$ |
| db38 | $${\color{green} \text{7.61e-17}}$$ | $${\color{green} \text{5.05e-16}}$$ | db4 | $${\color{green} \text{8.52e-17}}$$ | $${\color{green} \text{1.25e-16}}$$ | db5 | $${\color{green} \text{1.30e-16}}$$ | $${\color{green} \text{1.50e-16}}$$ |
| db6 | $${\color{green} \text{4.77e-17}}$$ | $${\color{green} \text{1.44e-16}}$$ | db7 | $${\color{green} \text{1.18e-16}}$$ | $${\color{green} \text{1.86e-16}}$$ | db8 | $${\color{green} \text{3.52e-17}}$$ | $${\color{green} \text{1.99e-16}}$$ |
| db9 | $${\color{green} \text{7.84e-17}}$$ | $${\color{green} \text{2.72e-16}}$$ | dmey | $${\color{red} \text{6.64e-01}}$$ | $${\color{red} \text{6.64e-01}}$$ | haar | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ |
| rbio1.1 | $${\color{green} \text{5.75e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ | rbio1.3 | $${\color{green} \text{5.80e-17}}$$ | $${\color{green} \text{7.05e-17}}$$ | rbio1.5 | $${\color{green} \text{5.84e-17}}$$ | $${\color{green} \text{7.48e-17}}$$ |
| rbio2.2 | $${\color{green} \text{3.11e-17}}$$ | $${\color{green} \text{2.08e-16}}$$ | rbio2.4 | $${\color{green} \text{3.90e-17}}$$ | $${\color{green} \text{6.03e-17}}$$ | rbio2.6 | $${\color{green} \text{9.92e-17}}$$ | $${\color{green} \text{1.22e-16}}$$ |
| rbio2.8 | $${\color{green} \text{1.49e-16}}$$ | $${\color{green} \text{1.09e-16}}$$ | rbio3.1 | $${\color{green} \text{9.97e-17}}$$ | $${\color{green} \text{2.27e-16}}$$ | rbio3.3 | $${\color{green} \text{8.47e-17}}$$ | $${\color{green} \text{1.75e-16}}$$ |
| rbio3.5 | $${\color{green} \text{2.20e-16}}$$ | $${\color{green} \text{1.31e-16}}$$ | rbio3.7 | $${\color{green} \text{1.91e-16}}$$ | $${\color{green} \text{1.44e-16}}$$ | rbio3.9 | $${\color{green} \text{2.48e-16}}$$ | $${\color{green} \text{2.58e-16}}$$ |
| rbio4.4 | $${\color{yellow} \text{2.74e-12}}$$ | $${\color{yellow} \text{2.74e-12}}$$ | rbio5.5 | $${\color{yellow} \text{1.09e-08}}$$ | $${\color{yellow} \text{1.09e-08}}$$ | rbio6.8 | $${\color{yellow} \text{3.81e-13}}$$ | $${\color{yellow} \text{3.81e-13}}$$ |
| sym10 | $${\color{yellow} \text{3.20e-14}}$$ | $${\color{yellow} \text{3.20e-14}}$$ | sym11 | $${\color{yellow} \text{1.13e-13}}$$ | $${\color{yellow} \text{1.13e-13}}$$ | sym12 | $${\color{yellow} \text{7.85e-14}}$$ | $${\color{yellow} \text{7.85e-14}}$$ |
| sym13 | $${\color{yellow} \text{5.74e-12}}$$ | $${\color{yellow} \text{5.74e-12}}$$ | sym14 | $${\color{yellow} \text{3.44e-12}}$$ | $${\color{yellow} \text{3.44e-12}}$$ | sym15 | $${\color{yellow} \text{2.83e-11}}$$ | $${\color{yellow} \text{2.83e-11}}$$ |
| sym16 | $${\color{yellow} \text{1.36e-11}}$$ | $${\color{yellow} \text{1.36e-11}}$$ | sym17 | $${\color{yellow} \text{6.70e-12}}$$ | $${\color{yellow} \text{6.70e-12}}$$ | sym18 | $${\color{yellow} \text{2.18e-10}}$$ | $${\color{yellow} \text{2.18e-10}}$$ |
| sym19 | $${\color{yellow} \text{4.16e-11}}$$ | $${\color{yellow} \text{4.16e-11}}$$ | sym2 | $${\color{yellow} \text{1.28e-12}}$$ | $${\color{yellow} \text{1.28e-12}}$$ | sym20 | $${\color{yellow} \text{2.03e-10}}$$ | $${\color{yellow} \text{2.03e-10}}$$ |
| sym3 | $${\color{yellow} \text{1.44e-10}}$$ | $${\color{yellow} \text{1.44e-10}}$$ | sym4 | $${\color{yellow} \text{1.64e-12}}$$ | $${\color{yellow} \text{1.64e-12}}$$ | sym5 | $${\color{yellow} \text{8.08e-13}}$$ | $${\color{yellow} \text{8.08e-13}}$$ |
| sym6 | $${\color{yellow} \text{2.39e-12}}$$ | $${\color{yellow} \text{2.39e-12}}$$ | sym7 | $${\color{yellow} \text{1.93e-12}}$$ | $${\color{yellow} \text{1.93e-12}}$$ | sym8 | $${\color{yellow} \text{1.56e-12}}$$ | $${\color{yellow} \text{1.56e-12}}$$ |
| sym9 | $${\color{yellow} \text{1.78e-14}}$$ | $${\color{yellow} \text{1.77e-14}}$$ |  |  |  |  |  |  |

The above table is colored as: `green` if $L_2 < 10^{-14}$, `yellow` if $L_2 < 10^{-7}$ and otherwise `red`.

![Error scatter](./figs/l2_scatter.svg)

We can see that FP64 error is highly correlated with Q-field error for most wavelet families, suggesting that coefficient conversion to `float64` is usually not the main source of error. The main error appears to arise earlier in the pipeline, likely during projection of the input matrix to determinant exactly $1$.

![Error scatter](./figs/l2_scatter_coif.svg)

Only for coiflets the datatype matters and converting fractions to `float64` increases error rate by a constant factor.

## 4. Comparing numerical stability


