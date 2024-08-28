[Meeting note](https://www.dropbox.com/s/ebpsm4r1vqf635s/Shih%20Huan%20Summer%20project.pdf?dl=0)

# [Protein aggregation](https://pubmed.ncbi.nlm.nih.gov/21842954/ "Cohen JCP")
### Master equation
$\frac{df(t,j)}{dt}=2m(t)k_+(f(t,j-1)-f(t,j))+2k_{off}(f(t,j+1)-f(t,j))-k_-(j-1)f(t,j))+2k_-\displaystyle{\sum_{i=j+1}^{\infty}}f(t,i)+k_2m(t)^{n2}M(t)\delta_{j,n_2}+k_nm(t)^{n_c}\delta_{j,n_c}$


### Equations for the moments (similar notation as Georg's thesis):
$\frac{dP(t)}{dt}=k_-[M(t)-(2n_c-1)P(t))]+k_2m(t)^{n_2}M(t)+k_nm(t)^{n_c}
\\\coloneqq \alpha_0P(t)+\beta_0(t)M(t)+c(t)$

$\frac{dM(t)}{dt}=2[m(t)k_+-k_{off}-k_-n_c(n_c-1)/2]P(t)+n_2k_2m(t)^{n_2}M(t)+n_ck_nm(t)^{n_c}
\\\coloneqq \alpha_1P(t)+\beta_1(t)M(t)+n_cc(t)$

### Solutions to the master equation
>*Assumption:* 
>1. *the monomer concentration is constant.*
>1. *monomer-dependent secondary nucleation is negligible*

$\frac{dP_0(t)}{dt}=k_-[M_0(t)-(2n_c-1)P_0(t))]+k_nm(0)^{n_c}$

$\frac{dM(t)}{dt}=2[m(0)k_+-k_{off}-k_-n_c(n_c-1)/2]P_0(t)+n_ck_nm(0)^{n_c}$

The solution is

$$P_0(t)=C_1e^{\kappa_1t}+C_2e^{\kappa_2t}-\frac{\eta_2}{\xi_2}$$

$$M_0(t)=\frac{C_1\xi_2}{\kappa_1}e^{\kappa_1t}+\frac{C_2\xi_2}{\kappa_2}e^{\kappa_2t}-\frac{\eta_1}{k_-}-\frac{\xi_1\eta_2}{\xi_2}$$
, where

$\kappa_{1,2}=\frac{1}{2}(-k_-\xi_1\pm\sqrt{k^2_-\xi^2_1+4k_-\xi_2})$

$\xi_1=2n_c-1$

$\xi_2=2m(0)k_+-2k_{off}-k_-n_c(n_c-1)$

$\eta_1=k_nm(0)^{n_c}$

$\eta_2=n_ck_nm(0)^{n_c}$

$C_{1,2}=\frac{1}{1-\frac{\kappa_{2,1}}{\kappa_{1,2}}}(\frac{\eta_2}{\xi_2}-\frac{\eta_1\kappa_{2,1}}{k_-\xi_2}-\frac{\xi_1\eta_2\kappa_{2,1}}{\xi_2^2}+P(0)-M(0)\frac{\kappa_{2,1}}{\xi_2})$

>*Here, we make a further assumption:*
>1. *$m(0)k_+>>k_-$*

Hence, 

$\xi_2\approx 2[m(0)k_+-k_{off}]$

$\kappa_{1,2}\approx \pm\frac{1}{2}\sqrt{k^2_-\xi^2_1+4k_-\xi_2}\approx \pm\sqrt{k_-\xi_2}=\pm\sqrt{2[m(0)k_+-k_{off}]k_-}\coloneqq\pm\kappa$

$P_0(t)=C_1e^{\kappa t}+C_2e^{-\kappa t}-\frac{n_ck_nm(0)^{n_c}}{2[m(0)k_+-k_{off}]}$

$M_0(t)=\frac{C_12[m(0)k_+-k_{off}]}{\kappa}e^{\kappa t}-\frac{C_22[m(0)k_+-k_{off}]}{\kappa}e^{-\kappa t}-\frac{k_nm(0)^{n_c}}{k_-}$

$C_{1,2}\approx \frac{1}{2}(P(0)\pm\frac{\kappa M(0)}{2[m(0)k_+-k_{off}]}\pm \frac{\kappa k_nm(0)^{n_c}}{2[m(0)k_+-k_{off}]k_-})$

### Analysis of limiting cases

### Analysis of the central moments

#### Mean filament length
$\mu(t)=\frac{M(t)}{P(t)}$

>At early times in *growth* phase, when the monomer concentration is held constant, denoted by a subscript 0,

$\mu_0(t)=\frac{2[k_+m(0)-k_off]}{\kappa}tanh(\frac{\kappa t}{2})$

For long time,

$\mu_0(\infty)=\frac{2[k_+m(0)-k_off]}{\kappa}=\frac{2[k_+m(0)-k_off]}{\sqrt{2[m(0)k_+-k_{off}]k_-}}=\sqrt{\frac{2[m(0)k_+-k_{off}]}{k_-}}$

>When $k_+m(0)>>k_{off}$

$\mu_0(\infty)=\sqrt{\frac{2m(0)k_+}{k_-}}$

> At later times, as the monomer is depleted, the full nonlinear solutions for M(t) and P(t ) show that the length decreases due to fragmentation dominating over [elongation](https://drive.google.com/file/d/1Zsb-pEJCJipKvLfSNm7JHBHA3KU5va7M/view?usp=sharing "time evolution of average length").


### Convert Gillespie to bulk rates
The units of concentation and rate constant in Gillespie are $_gC=(molecules)$ and $_gk=(time)^{-1}(molecule)^{-n}$. There is lack of volume term is unit. he relevant volume should be the cell volume $(10^{-5}m)^3\approx 1pl$.
#### Concentration
To convert Gillespie concentration to real concentration, we do the following calculation
$$C={_g}C*\frac{1}{VN_0}={_g}C*\frac{10^{12}}{6\times10^{23}}$$

According to experiments, concentrations of protein aggregates are in the range of $$10^{-9}M<C<10^{-5}M$$
Hence, $$600<{_gC}<6000000$$
#### Rate constant
Typical $k_+ =3\times10^{6}M^{-1}s^{-1}$, $k_+ = {_g}k*VN_0$. 

Hence, $_gk_+=5*10^{-6}s^{-1}(molec.)^{-1}$ 

$k_{off} =10^{-14}(sec)^{-1}$

Similarly, $_gk_n=5*10^{-16}s^{-1}(molec.)^{-1}$, $_gk_-=10^{-8}s^{-1}=k_-$, $k_2 = 10^4M^{-2}s^{-1}$ is equivalent to $k_- = 10^{-8}s^{-1}$ at $1\mu M$.

### Estimate the simulation time
From the [performance test](https://docs.google.com/spreadsheets/d/1CLMphbjoKtfBzVSIach01P74-mYVi8oHnucBW3zwJ0Y/edit?usp=sharing), to simulate a total fibril length of 20, it requires 100s.
To scale to a longer fibril, the time requires
$$c\times(20)^2=100s$$
$$c\times(10^5)^2=2\times10^9s\approx5.5\times10^5hr\approx63yr$$
for size = 10000 fibril. 

----

Simulation conditions:

$_gk_+=5*10^{-8}s^{-1}(molec.)^{-1}$

$_gk_n=5*10^{-16}s^{-1}(molec.)^{-1}$

$_gk_- = 10^{-6}s^{-1}$ at $1\mu M$.

$_gk_{off} =10^{-14}(sec)^{-1}$

$600<{_gC}=m_0<6000000$

Expected average length = $\sqrt{\frac{k_+m}{k_-}}= \sqrt{\frac{0.00000005*100000}{0.000001}}\approx 70$


$t = 10000000s, steps = 100$

Number of trajectories = $50$

Command line: 
> python Trial_of_fibril.py -kn 0.0000000000000005 -k+ 0.00000005 -k- 0.000001 -t 10000000 -s 100 -m *variable* -si *variable* -nc 2

Results:

|Test|k+|k-|kn|koff|Monomer|mu|simu_time|Time lapsed(s)|
|--|--|--|--|--|--|--|--|--|
|[1](https://drive.google.com/drive/folders/1tfLI7LaKOBCA7KZ02mXJ8B9CVynnTW5x?usp=sharing)|5e-07|1e-5|5e-12|1e-14|10000|31|10kappa^-1|82|
|[2](https://drive.google.com/drive/folders/1PSzZkhybYjH4cKAqt495gqDukM8HQreZ?usp=sharing)|5e-08|1e-6|5e-13|1e-14|10000|31|10kappa^-1|78|
|[3](https://drive.google.com/drive/folders/1wGQPBlnAReI3hkU3TecssuA7qowlsLqw?usp=sharing)|5e-08|1e-6|5e-15|1e-14|10000|31|100kappa^-1|73|
|[4](https://drive.google.com/drive/folders/1rH8IpQdwlV0gglCmh1GUn4TsqLpJzo4_?usp=sharing)|5e-08|1e-6|5e-16|1e-14|10000|31|100kappa^-1|43|
|[5](https://drive.google.com/drive/folders/1P3tyLJB9chWP5PCb9D1GzxJF3JsdrnwR?usp=sharing)|5e-08|1e-6|5e-16|1e-14|20000|44|10kappa^-1|129|
|[6](https://drive.google.com/drive/folders/12uOOCnkXwTl5msmmuuLSpDyvUGY1meXu?usp=sharing)|5e-08|1e-6|5e-16|1e-14|20000|44|100kappa^-1|176|

|Test|k+|k-|kn|koff|Monomer|mu|simu_time|Time lapsed(s)|
|--|--|--|--|--|--|--|--|--|
|[1](https://drive.google.com/file/d/1UvCaYXt1UvvNvfgn9x7DK_-Y8iwu0iqE/view?usp=sharing)|5e-07|1e-5|5e-12|1e-14|1000|10|30kappa^-1|4|
|[2](https://drive.google.com/file/d/12jvMo0A_zDzqzVCncRZhq9VVLwTPr7be/view?usp=sharing)|5e-07|1e-5|5e-12|1e-14|10000|31|30kappa^-1|82|
|[3](https://drive.google.com/file/d/1GqySia9hsL-M-e9Wd6NKg-eVddrltNCc/view?usp=sharing)|5e-07|1e-5|5e-12|1e-14|20000|44|30kappa^-1|117|


### Fibril length coarse-graining - bulk --> This part does not have a conclusion

### Fibril length coarse-graining - stochastic

