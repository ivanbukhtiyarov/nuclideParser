# OpenBSD
## Description
Nuclide kinetic solver
## Platforms

Platforms:

*   Linux
*   Mac OS X
*   Windows
## Requirenments
* a C++11-standard-compliant compiler
## External libraries

libraries:

*   Xtl vers 0.6.13 (https://github.com/xtensor-stack/xtl)
*   Xtensor vers 0.21.4 (https://github.com/xtensor-stack/xtensor)
*   Xtensor-blas vers 0.21.4 (https://github.com/xtensor-stack/xtensor-blas)
*   Google-test vers 1.10.x (https://github.com/google/googletest)
## Getting started
### Config
Config default
``` xml
<?xml version="1.0"?>
<configure>

  <chain>./Xmls/chain_simple.xml</chain>
  <reaction>./Xmls/reactions.xml</reaction>
  <numbers>1</numbers>
  <timestep>720000</timestep>

</configure>
```
tag **chain** : path to xml file with nuclide chain (examples in Xmls dir `chain_simple.xml`)

tag **reaction** : path to xml file with nuclide reactions(examples in Xmls dir `reactions.xml`)


## Theory
   
The system of equations nuclide kinetics is inhomogeneous linear system of equations:


$$
\frac{d y_{i}}{d t}=\sum_{j=1}^{n} a_{i j} y_{j}+q_{i}, \quad(i=1,2, \dots, n)     \quad  \quad \quad   (1)
$$

with the initial condition:

$$
y_{i}(0)=y_{i 0}, \quad(i=1,2, \dots, n) \quad  \quad \quad \quad \quad\quad \quad (2)
$$

where y<sub>i</sub> is the concentration of the i th nuclide, a<sub>ij</sub>  are the coefficients characterizing the channels of the transformations of the i th nuclide from the j–th nuclide, qi , an external source, n is the number of nuclides, y<sub>i0</sub>  is the concentration of the i th nuclide at time t<sub>0</sub> .

- ## Iteration Method 
  
The system of equations (1) – (2) in matrix form:
$$
\frac{d \vec{y}}{d t}=\hat{A} \vec{y}+\vec{q}, \quad\left(A \equiv\left[a_{i j}\right]\right)\quad \quad \quad\quad \quad\quad \quad  (3)
$$
$$
\vec{y}(0)=\vec{y}_{0} \quad \quad \quad\quad \quad\quad \quad \quad \quad \quad\quad \quad\quad \quad \quad (4)
$$
The system of equations (1) – (2) can be solved analytically as the sum of its particular solution and the General solution of the corresponding homogeneous linear system of equations. In this case, there is a need matrix of size n × n multiplication. Therefore, this method of solution is often used for a partial transition matrices, for example, triangular [18, 21].
Another widely known method of solving the problem is to decompose the solution in the range of the exponential function, which leads to the need to use a recurrence relation [11].
In the code the solution of equations (1) – (2) is the iterative method. Detailed description of the method can be found in the papers [22 – 24]. For the k th iteration solution is obtained:
$$
y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)+\sum_{j \neq i} d y_{k-1}^{j} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}, \quad(j=1,2, \ldots, k) ,(5)
$$
$$
d y_{k}^{i}=\sum_{j \neq i} d y_{k-1}^{j} \frac{\lambda^{j-i}}{\lambda_{p}^{j}}\left[1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right], \quad(j=1,2, \ldots, k),(6)
$$
where τ is the time step; λ<sup>j→i</sup> is the rate of formation of the i th nuclide from the j th nuclide, taking into account the probability of such a process (with the possibility of branching); and the speed of withdrawal of the nuclide j and i respectively at the expense of all processes.

The final solution is:

$$
\begin{aligned}
&y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)+\\
&\begin{array}{l}
+\sum_{j_{k} \neq i} \frac{\lambda^{j_{k} \rightarrow i}}{\lambda_{p}^{j_{k}}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)^{k-1} \prod_{m=1}^{k-1}\left(\sum_{j_{n} \neq i} \frac{\lambda^{j_{m} \rightarrow i}}{\lambda_{p}^{j_{n}}}\right) d y_{0}^{j_{k}} \approx \\
\approx y_{0}^{i}(\tau)+\sum_{j \neq i}\left\{d y_{0}^{j} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\left[1+\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right) \sum_{j \neq i} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}}+\right.\right.
\end{array}\\
&\left.\left.+\ldots+\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)^{k-1} \prod_{m=1}^{k-1}\left(\sum_{j_{n} \neq i} \frac{\lambda^{j_{n}-i}}{\lambda_{p}^{j_{n}}}\right)\right]\right\}  (7)
\end{aligned}
$$
Equation (5) describes the concentration of a nuclide, given the possible departure newcomers nuclei of that nuclide at the end of the time step. Equation (6) takes into account the formation of new nuclei of a nuclide of the other, broken at the k-th iteration.
The algorithm takes into account two main types of channels of nuclear transformations: radioactive decay of nuclei and nuclear reactions initiated by neutrons.

The rate of nuclear reactions for the i-th nuclide is described by the expression:

$$
\lambda_{r}=\frac{1}{V} \iint_{V_{F}} \sigma_{r}(E, \vec{r}) \varphi(E, \vec{r}) d E d V, (8)
$$

where λ is the reaction rate, V – volume of the computational cell, σ is micromachine reactions caused by neutrons with energy E in the point φ is the flux density of neutrons.
For fission expression (8) takes the form:
$$
\lambda_{f i}=\frac{1}{V} \iint_{V E} \sigma_{f i}(E, \vec{r}) \varphi(E, \vec{r}) d E d V ,(9)
$$
For speed capture process λ<sub>ci</sub>:
$$
\lambda_{c i}=\frac{1}{V} \iint_{V E} \sigma_{c i}(E, \vec{r}) \varphi(E, \vec{r}) d E d V ,(10)
$$
When calculating the speeds of neutron reactions assumes the immutability of the absolute density of the neutron flow on the time interval, i.e. throughout step τ used constant speed processes.
The transfer speed of the i th nuclide in the j th nuclide in the radioactive decay with half-lives T <sub>½</sub> and the probability of decay channel of the ε<sub>i→j</sub> is described by the expression:
$$
\lambda^{i \rightarrow j}=\varepsilon_{i \rightarrow j} \frac{\ln (2)}{T_{1 / 2}^{i}} ,(11)
$$
Given the properties of the expression
$$
\begin{aligned}
&\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)<1, \quad \lambda_{p}^{i} \rightarrow \infty\\
&\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right) \rightarrow 0, \quad \lambda_{p}^{i} \rightarrow 0
\end{aligned} (12)
$$
you can verify that the formation of "new" nuclei in (6) at the end of the time step tend to zero, which ensures convergence of the iteration process.
Moreover, the solution of (5) – (6) is always non-negative, i.e. the solution of the system of equations (1) – (2) y<sub>i</sub> exists and it is always positive, since in the sum (5) is always at least one summand is different from zero because the initial concentration of at least one nuclide is positive (otherwise, the task loses physical meaning).
The iterative process of solving equations (5) – (6) continues until the condition is met:
$$
\left|1-\frac{y_{n}^{k}}{y_{n-1}^{k}}\right| \leq \delta_{\max }, \quad \forall k  \quad  \quad  \quad  \quad  \quad  \quad  \quad  (13)
$$
where δ<sub>max</sub> is the maximum allowable value of the variation of the concentrations of nuclides at two adjacent iterations specified by the user in the calculation options. The accuracy of the iterative solution the default is 10-5, but can be changed by the user.
From equation (7) it follows that when k → ∞ the last term in the series is actually a product of two multiplicands, each of which is raised to the power of (k - 1) by the number of obviously smaller units and tend to zero, which ensures convergence. According to [25] a number of (7) can be interpreted as the Neumann series, convergence has been proven.
During the iterative process operated with only non-negative values of equations (5) and (6). Thus, a solution exists, it is positive, because the equation (5) is always one summand is different from zero.
The calculation of the decay heat WOST is inextricably linked with the concentration of radioactive nuclides in the material and the energy released during radioactive decay of the nucleus, and is determined by the formula:
$$
W_{o c m}=1,60219 \cdot 10^{-13} \cdot \sum_{j} y_{j} \lambda_{j} E_{j} ,\quad \quad (14)
$$
where Ej is the heat generation by the decay of nuclide j, MeV/dis.; λj – decay constant of j-th radioactive nuclide with-1; 1,60219 ∙ 10-13 – the conversion factor from MeV to watts.
The calculation of activity of nuclide Aj(t) at time t is a product of its nuclear density yi(t) by a constant decay in the whole volume of the computational cell:
$$
A_{j}(t)=\lambda_{j} \int_{0}^{V} y_{j}(\vec{r}) d v,\quad \quad \quad \quad \quad (15)
$$
where r - the coordinate of the calculated point, V - its volume.

- ## Integration Method 
  
The nonau- tonomous linear ODE appears a number of times as a subproblem that needs to be solved. The Magnus expansion is an infinite series solution to. The Magnus expansion:
$$
\begin{aligned}
y(t) &=\exp \Omega(t) y(0) \\
\Omega(0) &=\mathbf{O} \\
\Omega(t) &=\sum_{k=1}^{\infty} \Omega_{k}(t)
\end{aligned} (16)
$$
Essentially, the solution to 𝑦(𝑡) is given by a single matrix exponential of the matrix Ω(𝑡). This matrix is formed via an infinite sum of terms, of which the first three are given below. Here, [·, ·] is the matrix commutator, where [𝐴, 𝐵] = 𝐴𝐵 − 𝐵𝐴.
$$
\begin{array}{l}
\Omega_{1}(t)=\int_{0}^{t} A\left(t_{1}\right) d t_{1} \\
\Omega_{2}(t)=\frac{1}{2} \int_{0}^{t} \int_{0}^{t_{1}}\left[A\left(t_{1}\right), A\left(t_{2}\right)\right] d t_{2} d t_{1} \\
\Omega_{3}(t)=\frac{1}{6} \int_{0}^{t} \int_{0}^{t_{1}} \int_{0}^{t_{2}}\left(\left[A\left(t_{1}\right),\left[A\left(t_{2}\right), A\left(t_{3}\right)\right]\right]+\left[A\left(t_{3}\right),\left[A\left(t_{2}\right), A\left(t_{1}\right)\right]\right]\right) d t_{3} d t_{2} d t_{1}
\end{array} (17)
$$
Further terms become increasingly complex. One way to transform this into an effective integrator is to truncate the sum in Equation (16). If accuracy is still too low, one could perform substepping in which, instead of integrating in one shot to 𝑡 + h, 𝑚 integrations of h/𝑚 are performed. If one performs this substepping with Ω truncated to Ω1, the result is Equation (18). This happens to be identical to the substepping method used for the CE/LI and LE/QI predictor-corrector methods.
$$
\begin{aligned}
A_{s} &=\int_{t_{n}+\frac{s-1}{m}}^{t_{n}+\frac{s}{m} h} A(s) d s \\
y\left(t_{n}+h\right) &=\operatorname{expm}\left(A_{m}\right) \operatorname{expm}\left(A_{m-1}\right) \ldots \operatorname{expm}\left(A_{1}\right) y\left(t_{n}\right)
\end{aligned}(18)
$$ 
A much cheaper approach than directly evaluating the expansion is to form a quadrature. There are a wide variety of quadrature options. One option, Equation (19), only requires one matrix exponential, but requires matrix commutators in order to do.
$$
\begin{aligned}
c &=\frac{1}{2} \mp \frac{\sqrt{3}}{6} \\
A_{i} &=A\left(t+c_{i} h\right) \\
\Omega^{[4]}(h) &=\frac{h}{2}\left(A_{1}+A_{2}\right)-\frac{h^{2} \sqrt{3}}{12}\left[A_{1}, A_{2}\right] \\
y(t+h) &=\exp \left(\Omega^{[4]}(h)\right) y(t)
\end{aligned}(19)
$$
Unfortunately, this particular method proved to be unstable during testing. An alternative form, and one that did not have such an issue, is the commutator free integrator shown in Equation (20). This form removes the need to compute commutators in exchange for requiring multiple matrix exponentials. This algorithm will be abbreviated as “CFQ4” in the rest of the paper.
$$
\begin{aligned}
c &=\frac{1}{2} \mp \frac{\sqrt{3}}{6} \\
a &=\frac{1}{4} \pm \frac{\sqrt{3}}{6} \\
A_{i} &=h A\left(t+c_{i} h\right) \\
y(t+h) &=\exp \left(a_{2} A_{1}+a_{1} A_{2}\right) \exp \left(a_{1} A_{1}+a_{2} A_{2}\right) y(t)
\end{aligned}(20)
$$
It is unclear why one works and the other does not, but there are at least a few cases where the use of commutators reduces the stability of numerical integration. Additionally, methods based on the Magnus expansion directly may fail if the expansion does not converge. This can happen if
$$
\int_{0}^{t}\|A(s)\|_{2} d s \geq \pi
$$
The other problem is that the eigenvalues of Ω are not known very rigorously. This can cause issues with the Chebyshev rational approximation matrix exponential recommended for use with depletion.

- ## Chebyshev Rational Approximation 
The Chebyshev rational approximation method (CRAM) is a relatively straight- forward algorithm. A rational function 𝑟^𝑘,𝑘(𝑥) is found that minimizes the max- imum error with regard to the scalar exponent along the negative real axis
The defining equation is Equation (21), where 𝜋𝑘,𝑘 is the set of all rational functions with numerators and denominators of order 𝑘. As 𝑘 increases, the accuracy of the approximation also increases.
$$
\sup _{x \in \mathbb{R}_{-}}\left|\hat{r}_{k, k}(x)-e^{x}\right|=\inf _{r_{k, k} \in \pi_{k, k}}\left\{\sup _{x \in \mathbb{R}_{-}}\left|r_{k, k}(x)-e^{x}\right|\right\} (21)
$$
Once the function 𝑟^𝑘,𝑘(𝑥) is known, it can be rearranged to reduce costs further or to improve numerical stability. The incomplete partial fraction (IPF) form, shown in Equation (22), is a good combination of numerical stability and efficiency. The values 𝛼𝑙 and 𝜃𝑙 are tabulated and are available for a variety of values of 𝑘 up to 48. In the IPF form, only sparse matrix solves are necessary to compute the action on a vector.
$$
\hat{r}_{k, k}(x)=\alpha_{0} \prod_{l=1}^{k / 2}\left(1+2 \Re\left\{\frac{\tilde{\alpha}_{l}}{x-\theta_{l}}\right\}\right) (22)
$$
CRAM is both efficient and highly accurate over the domain in which it is derived. However, eigenvalues with extremely large imaginary components or positive real components will reduce the accuracy. As such, CRAM is not recommended for use in highly oscillatory problems or those with possible exponential growth such as reactor dynamics.
