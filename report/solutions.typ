#import "@preview/cetz:0.3.4"
#import "@preview/fletcher:0.5.7" as fletcher: diagram, node, edge, shapes
#import "@local/uva-game-theory:0.1.0": *

#set math.equation(numbering: "(1)", number-align: bottom)
#set heading(numbering: "1.1")
#set enum(numbering: "Q1.i.")
#set block(width: 100%)

#show heading.where(level: 3): set heading(hanging-indent: 0em, numbering: (..nums) => {
  let numbering = {
    if nums.at(1) == 0 {
      nums.at(2)
    } else {
      nums.pos().slice(1).map(str).join(".")
    }
  }
  [*Q#{numbering}.*]
})
#show heading.where(level: 3): it => {
  set text(weight: "regular")
  set par(first-line-indent: 4em)
  it
}

#let microstate(x) = $bold(#x)$
#let ensembleavg(varname, microstatevar: $s$, probfn: $P$) = $sum_(microstate(microstatevar)) varname(microstate(microstatevar)) dot probfn(microstate(microstatevar))$

= Modelling the activity of a single neuron

#let neuron_results = json("../results/neuron_modelling.json")

=== Can you plot the distribution $P(tau)$ of the time intervals $tau$ between successive spikes? Check that there is indeed a rafactory period, i.e., a time interval $tau_0$ after each spike, during which the neuron doesn't spike again. What is the duration $tau_0$ for this time interval?
  #solution[
  @fig:neuron-empirical-interspike-distribution shows the distribution of inter-spike durations in the dataset as a histogram with 50 bins of uniform width. We obseve a refractory period of $tau_0 = #calc.round(neuron_results.tau0, digits: 2)"ms"$, indicated by the lack of inter-spike durations at the left-hand side of the distribution, in the interval $[0, tau_0)$.
  #figure(
    image("../results/figures/neuron_interspiking_distribution.svg"),
    caption: [Distribution over the duration between activity spikes for a neuron. A refractory period of $tau_0 approx 1.9 m s$ is identified as the minimum observed duration between spikes.]
  ) <fig:neuron-empirical-interspike-distribution>
  ]

=== Can you check that the decay of the distribution $P(tau)$ of inter-spike intervals is indeed exponential? Measure the corresponding decay rate $lambda$

#solution[
  We say $P(tau)$ decreases with exponential decay for increasing $tau > tau_0$ if it is well-described by a model of the form
  $
  P(tau - tau_0) = P(tau_0) e^(-lambda (tau - tau_0))
  $ <eq:exponential-decay>
  To test this hypothesis, we can check whether the observed inter-spike durations display a good linear fit under the equivalent linear model:

  $
  P(tau - tau_0) &= &&P(tau_0) e^(-lambda (tau - tau_0)) \
  log P(tau - tau_0) &= &&log P(tau_0) - lambda (tau - tau_0) \
  -log P(tau - tau_0) &= -&&log P(tau_0) + lambda (tau - tau_0)
  $<eq:exponential-decay-linear-form>

  We estimate $P(tau - tau_0)$ by subtracting the calculated refractory period from each observed inter-spike duration, then binning the data into 50 bins of uniform width. We normalise the observation count in each bin with respect to the total number of observations ($n=#neuron_results.n$), and take this to be the empirical probability for the value of $tau$ at the center of each bin. 

  @fig:neuron-log-linear-fit shows the empirical probabilities transformed using @eq:exponential-decay-linear-form. We observe a good linear fit, indicating exponential decay in $P(tau)$ with $lambda = #calc.round(neuron_results.exponential_fit.lambda, digits: 4)$ given by the slope of the linear fit. The quality of the linear fit is reduced for large $tau$; however, this is explained by a smaller sample size in this range.

  #figure(
    image("../results/figures/neuron_exponential_fit.svg"),
    caption: [Neuron inter-spike durations (without refractory period) shows a linear relationship under a log-transform, with slope $lambda = #calc.round(neuron_results.exponential_fit.lambda, digits: 4)$ and intercept $P(tau_0) = #calc.round(-calc.ln(neuron_results.exponential_fit.p_tau_0), digits: 4)$. $R=#calc.round(neuron_results.exponential_fit.r_value, digits: 4)$.]
  ) <fig:neuron-log-linear-fit>
]

=== Can you deduce an analytical expression for the distribution of inter-spike time interval $P(tau)$ of the delayed Poisson process as a function of $lambda$ and $tau_0$? Compare your model distribution to the one obtained from the data.

#solution[
  In the previous question we found that $P(tau)$ exhibited exponential decay for $tau > tau_0$. As such, we can model $P(tau)$ as a shifted exponential distribution where $P(tau) = 0$ for $tau < tau_0$, and the rate $lambda$ is as calculated in *Q2*.
  $
  P(tau) = cases(
    0 quad &tau < tau_0 \
    lambda exp(-lambda (tau - tau_0)) quad &tau >= tau_0
  )
  $ <eq:analytical-interspike-duration>

  For @eq:analytical-interspike-duration to be a probability distribution it remains to show that it is well-normalised. To see that this is true, observe that the probability density for $tau < tau_0$ is 0, and the density for $tau >= tau_0$ is 1 (since this is just a standard exponential distribution). Thus it follows that 
  $ 
  integral_0^infinity P(t) dif t &= cancel(integral_0^(tau_0) P(t) dif t) + integral_(tau_0)^infinity lambda e^(-lambda (t - tau_0)) dif t  \
  &= integral_0^infinity lambda e^(-lambda t) dif t \
  &= 1
  $

  And so @eq:analytical-interspike-duration is an analytical expression for the distribution $P(tau)$.
]

=== Using your model, can you generate another 1000 (spike times) datapoints?

#solution[
  To sample 1000 new spike times, it suffices to sample 1000 inter-spike durations $bold(tau) = (tau_(n+1), ..., tau_(n+1000))$ from @eq:analytical-interspike-duration, and take the spike times $bold(t) = (t_(n+1), ..., t_(n+1000))$ with $t_(n+k) = t_n + sum_(i=1)^k tau_(n+i)$ (i.e., cumulative sum of interspike durations to $t_(n+k)$).

  However, the task of sampling from @eq:analytical-interspike-duration remains. We do this via the inverse-transform method. Letting $F(T)$ be the cumulative distribution associated with $P(tau)$, we sample $u in cal(U)_([0,1])$ and calculate $tau = F^(-1)(u)$ (so long as $F$ is invertible). Then each $tau$ is sampled with probability equal to $P(tau)$. $F(T)$ is derived as follows:
  $
  F(T) &= integral_0^T P(t) dif t \
  &= integral_(tau_0)^T P(t) dif t \
  &= lr([-exp(-lambda (t - tau_0))]_(tau_0)^T) \
  &= cases(
    0 &"if" T < tau_0 \
    1 - e^(-lambda (T - tau_0)) wide &"if" T >= tau_0
  )
  $ <eq:cumulative-distribution-interspike-duration>

  While $F(T)$ is piecewise this is not a problem for the inverse-transform since (in practice) $u != 0$. Indeed @eq:cumulative-distribution-interspike-duration is invertible for $u in (0,1)$, with:

  $
  u &= 1 - e^(-lambda (t - tau_0)) \
  1 - u &= e^(-lambda (t - tau_0)) \
  -lambda (t - tau_0) &= ln(1 - u) \
  t &= -1/lambda ln(1 - u) + tau_0
  $ <eq:inverse-transform-interspike-duration>

  @fig:additional-interspike-duration-samples shows the distribution over additional samples $bold(tau) = (tau_(n+1), ..., tau_(n+1000))$ as sampled using the inverse-transform method with @eq:inverse-transform-interspike-duration, overlaid on the distribution of observed data. 

  #figure(image("../results/figures/neuron_compare_empirical_theoretical.svg")) <fig:additional-interspike-duration-samples>
]

=== What is the average spiking rate $f$ of the neuron in the data? How is $f$ analytically related to $tau_0$ and $lambda$ that you have previously measured?

#solution[
  The average spiking rate $f$ is defined as the inverse of the expected duration between spikes, $E[tau]$. To determine the analytical form of $f$, we first derive $E[tau]$:
  $
  E[tau] &= cancel(integral_0^(tau_0) t P(t) dif t) + integral_(tau_0)^infinity t dot lambda e^(-lambda (t - tau_0)) dif t \
  &= [-t e^(-lambda (t - tau_0))]_(tau_0)^infinity - integral_(tau_0)^infinity -e^(-lambda (t - tau_0)) dif t wide "(integration by parts)"\
  &= lr([-t e^(-lambda (t - tau_0)) - (1/lambda) e^(-lambda (t - tau_0))])_(tau_0)^infinity \
  &= tau_0 + 1/lambda
  $ <eq:expected-interspike-duration>

  We then derive $f$ as $1/E[tau]$:
  $
  f = 1\/E[tau] &= 1/(tau_0 + 1/lambda) quad ==> quad lambda/(lambda tau_0 + 1)
  $ <eq:expected-spiking-rate>

  Evaluating @eq:expected-spiking-rate using our estimates for $tau_0$ (*Q1.*) and $lambda$ (*Q2.*), we obtain an estimate for the average spiking rate $f = #{calc.round(1000 * neuron_results.exponential_fit.lambda/(neuron_results.exponential_fit.lambda * neuron_results.tau0 + 1), digits: 2)}"hz"$. We can compare this to the average spiking rate as calculated directly from that data as the average value of $tau^(-1)$, $f_"data" = #calc.round(1000 * neuron_results.mean_spike_rate, digits: 2)"hz"$. The difference between the two values is likely the result of numerical error differences between the two methods.
]

= Modelling binary data with the Ising model

== Pairwise spin model
It is common to model the collective behavior of systems of binary variables with Ising-like models. To do so, we assume that the system is in a stationary state, and therefore that the datapoints are independently sampled from the same stationary probability distribution. We take this probability distribution to have the general form of an Ising model:

$
p_bold(g) (bold(s)) = 1/Z(bold(g)) exp(sum_(i=1)^n h_i s_i + sum_("pair"(i,j)) J_(i j) s_i s_j)
$ <eq:pairwise-ising-model>

where $n$ is the number of spin variables, $"pair"(i,j)$ denotes a summation over all possible pairs of distinct spin variables, #box[$bold(g) = (h_1, ..., h_n, J_(1,2), ..., J_(n-1,n))$] is a vector of (real) parameters, and $Z(bold(g))$ is a normalisation factor. There are several differences compared to the Ising model we have seen in class:

- There is a different external field $h_i$ for each spin $s_i$, which can take any real value. In particular, the $h_i$’s are not necessarily all positive or all negative.

- There is a different coupling parameter $J_(i j)$ for each pair of spins $s_i$ and $s_j$. The parameter $J_(i j)$ parametrises the strength of the coupling between $s_i$ and $s_j$, and can take any real value. In particular, the $J_(i j)$’s are not necessarily all positive or all negative.

- The summation is over all possible pairs ($i, j$) of spins, and not just over the “nearest neighbors”. The reason is that, in a general dataset, we have a priori no idea if there exists an underlying structure between the variables and if so, what that structure is, and therefore we don’t know which variables are “nearest neighbors”.

The general goal of the problem is to infer the set of parameters $bold(g) = (h_1, ..., h_n, J_(1,2), ... J_(n-1, n))$ that is the most appropriate to model the data, i.e. to find the parameters $bold(g)$ for which the probability distribution in @eq:pairwise-ising-model best fits the data. This way, we would infer from the data the underlying structure between the variables, and find which variables tend to be strongly influenced by an external parameter, and which pairs of variables tend to be more strongly coupled.

=== How many terms are in the sum over the $"pair"(i, j)$? Can you deduce what is the number of parameters in the vector $g = (h_1,..., h_n, J_(1,2),..., J_n)$? Can you re-write the sum over the $"pair"(i, j)$ as a double sum over $i$ and $j$ (without counting twice each pair)?
#solution[
Assuming that $"pair"(i,j)$ is order-independent, the number of pairs in the sum over $"pair"(i, j)$ is given by the number of ways which we can choose two distinct spins from a set of $n$ total spins, such that the order of the chosen spins doesn't matter: $n(n - 1)\/2$. So the total number of parameters in $bold(g)$ is:

$
n + (n (n-1))/2 = (n(n+1))/2
$

We can rewrite the sum using a double summation over $i$ and $j$ by enumerating the $k = n-1$ ways to choose the first spin, and then the $n - k$ ways to choose the second (given the first):

$
sum_("pair"(i,j)) J_(i,j) s_i s_j = sum_(i=1)^(n-1) sum_(j=i+1)^n J_(i,j) s_i s_j
$
]

=== Can you write down explicitly the terms in the exponential of @eq:pairwise-ising-model for a system with $n = 3$ spins?
#solution[
$
h_1 s_1 + h_2 s_2 + h_3 s_3 + J_(1,2) s_1 s_2 + J_(1,3) s_1 s_3 + J_(2,3) s_2 s_3
$
]

=== In @eq:pairwise-ising-model, we can recognize the Boltzmann distribution, in which the parameter $beta = 1\/(k_b T)$ was taken equal to 1 (more precisely, the constant $k_B$ was taken equal to 1, and the temperature parameter $T$ was absorbed in the parameters $h_i$ and $J_(i j)$). What is the energy function associated with the Boltzmann distribution in that case? What is the partition function and what is its general expression?
#solution[
The Boltzmann distribution has the form $P(hat(bold(s))) = 1/Z exp(-beta E(hat(bold(s))))$, where $beta = 1/(k_b T)$. Then from @eq:pairwise-ising-model we have
$
-1/(k_b T) E(hat(bold(s))) &= sum_(i=1)^n h'_i s_i + sum_("pair"(i,j)) J'_(i,j) s_i s_j \
==> quad E(hat(bold(s))) &= -k_b T lr([sum_(i=1)^n h'_i s_i + sum_("pair"(i,j)) J'_(i,j) s_i s_j]) \
&= -sum_(i=1)^n (T dot h'_i) s_i - sum_("pair"(i,j)) (T dot J'_(i,j) ) s_i s_j
$

As stated, the temperature parameter $T$ is absorbed into the parameters $h_i$ and $J_(i j)$ -- that is, for constant temperature $T$, we redefine $h_i = T dot h'_i$ and $J_(i j) = T dot J'_(i j)$, giving us the following energy function:

$
  E(hat(bold(s))) = -sum_(i=1)^n h_i s_i - sum_("pair"(i,j)) J_(i j) s_i s_j
$ <eq:ising-energy-function>

The partition function $Z$ is the normalisation factor for the Boltzmann distribution:

$
&sum_(hat(bold(s))) 1/Z(bold(g)) exp(-beta E(hat(bold(s)))) = 1 wide ==> wide Z(bold(g)) = sum_(hat(bold(s))) exp(-beta E(hat(bold(s))))
$

For the energy function in @eq:ising-energy-function $Z(bold(g))$ has the form:

$
Z(bold(g)) = sum_(hat(bold(s))) exp(sum_(i=1)^n h_i s_i + sum_("pair"(i,j)) J_(i j) s_i s_j)
$
]

=== Take a spin $s_i$: if $h_i$ is positive, which direction will $s_i$ tend to turn to, i.e., which direction of $s_i$ will minimize the associated energy $-h_i s_i$? Take a pair of spins $s_i$ and $s_j$: if $J_(i j)$ is positive, which configurations of $(s_i, s_j)$ minimize the coupling energy $-J_(i j) s_i s_j$?

  Assume that we have inferred the best parameters $h_i$ and $J_(i j)$ for the US supreme court dataset discussed in section 2. How would you interpret the sign of the inferred parameters $h_i$ and $J_(i j)$ in this context?

#solution[
Consider a spin $s_i$. If $h_i$ is positive, then the component of the spin's energy attributable to $h_i$ is minimised when #box[$s_i > 0$], such that $-h_i s_i < 0$. Similarly, take two spins $s_i, s_j$ with $i != j$, then if $J_(i j) > 0$ the energy attributable to the spins' interaction is minimised when $"sign"(s_i) = "sign"(s_j)$, such that $s_i s_j > 0$ and $-J_(i j) s_i s_j < 0$. 

If $h_i$ or $J_(i j)$ were negative then the opposite result holds ($s_i < 0$, or $"sign"(s_i) != "sign"(s_j)$ respectively). If $h_i = 0$ then $s_i$ has no preferred direction, i.e., the energy due to the field is minimised for any $s_i$. Likewise, if $J_(i j) = 0$ then any configuration of $s_i$ and $s_j$ minimises their interaction energy.

Suppose now that we have inferred the optimal parameters $h_i$ and $J_(i j)$ for the US supreme court dataset. We can interpret the sign of each parameter by considering its effect _in absence of the other's effect_. For instance, $"sign"(h_i)$ signifies the *political leaning* of the $i$'th judge's votes, in the absence of interactions with other judges. Analogously, $"sign"(J_(i j))$ signifies the *tendency for $i$ and $j$ to vote identically*, in the absence of their individual political leanings. 

Since we take $+1$ to represent a conservative vote and $-1$ a liberal one, $h_i < 0$ indicates that $i$ has a tendency to vote liberally, and vice versa. Judges $i$ and $j$ tend to vote similarly if $J_(i j) > 0$, and differently if $J_(i j) < 0$. 

Finally, $h_i = 0$ indicates no particular tendency to vote conservative or liberal, and $J_(i j) = 0$ implies no correlation between $i$ and $j$, aside from that which may arise from nonzero $h_i, h_j$.
]

== Observables

The important observables of the system are the local average magnetisations $angle.l s_i angle.r$ and the local correlations #box[$c_(i j) = angle.l s_i s_j angle.r - angle.l s_i angle.r angle.l s_j angle.r$], where the angle brackets $angle.l A(bold(s)) angle.r$ denote the ensemble average (or thermal average) of the microscopic quantity $A(bold(s))$.

=== Given a stationary probability distribution of the state $p_g (bold(s))$, what are the definitions of $angle.l s_i angle.r$ and of $angle.l s_i s_j angle.r$?
#solution[
For clarity define $s_i$ as a function which extracts the $i$'th element of a vector $s$, that is $s_i: bold(s) mapsto (bold(s))_i$. Then $angle.l s_i angle.r$ is the average local magnetisation of the $i$'th spin:

$
angle.l s_i angle.r = ensembleavg(s_i, probfn: p_bold(g))
$ <eq:avg-local-magnetisation>

And $angle.l s_i s_j angle.r$ is the average local correlation between the $i$'th and $j$'th spins:

$
angle.l s_i s_j angle.r = sum_microstate(s) s_i (microstate(s)) s_j (microstate(s)) dot p_bold(g)(microstate(s))
$ <eq:avg-local-correlation>
]

=== Consider a dataset $hat(bold(s))$ composed of $N$ independent observations of the spins: $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$. Let us denote by $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ the empirical averages of $s_i$ and of $s_i s_j$ respectively (i.e., their average values in the dataset). How would you compute $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ from the data?
#solution[
  Let $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$ be a dataset of observed microstates. We formalise the notion of the data distribution $p_hat(bold(s))$ as the proportion of $hat(bold(s))$ comprising observations of a particular microstate:

  //Let $D$ denote the 'dataset', i.e., the set of observation ${microstate(s)^((1)), ..., microstate(s)^((N))}$. We formalise the notion of the data distribution $P_D$ as the proportion of $D$ comprising observations of a particular microstate:
  
  $
  p_hat(bold(s)) (microstate(s)) = (\# {microstate(s)^((k)) in hat(bold(s)) | microstate(s)^((k)) = microstate(s)})/N
  $ <eq:data-probability>

  We define $angle.l s_i angle.r_D$ using $p_hat(bold(s))$, in such a form that it can be computed from $hat(bold(s))$:

  $
  angle.l s_i angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) dot p_hat(bold(s)) (microstate(s)) \
  &= cancel(sum_(microstate(s) in.not hat(bold(s))) s_i (microstate(s)) dot 0) + sum_(microstate(s) in hat(bold(s))) s_i (microstate(s)) dot p_hat(bold(s)) (microstate(s)) \
  &= 1/N sum_(microstate(s) in hat(bold(s))) s_i (microstate(s)) dot \# {microstate(s)^((k)) in hat(bold(s)) | microstate(s)^((k)) = microstate(s)} \
  &= 1/N sum_(k=1)^N s_i (microstate(s)^((k)))
  $ <eq:local-magnetisation-using-pdata>


  i.e., The empirical average local magnetisation of spin $i$ can be calculated as the average value of $s_i$ as observed in the data. Note that the sum index, $microstate(s) in hat(bold(s))$ on the second and third lines is taken to mean "$microstate(s)$ occurs in the dataset", rather than being an enumeration of the rows. 

  #block(breakable: false)[
  Likewise for $angle.l s_i s_j angle.r_D$:
  $
  angle.l s_i s_j angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) s_j (microstate(s)) dot p_hat(bold(s)) (microstate(s)) \
  &= 1/N sum_(microstate(s) in hat(bold(s))) s_i (microstate(s)) s_j (microstate(s)) dot \# {microstate(s)^((k)) in hat(bold(s)) | microstate(s)^((k)) = microstate(s)} \
  &= 1/N sum_(k=1)^N s_i (microstate(s)^((k))) s_j (microstate(s)^((k)))
  $ <eq:local-correlation-using-pdata>
  ]
]

=== Assume that the data is stationary and that each datapoint has been randomly sampled from $p(bold(s))$. Can you show that the empirical averages, $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$, converge to the model averages, respectively $angle.l s_i angle.r$ and $angle.l s_i s_j angle.r$, as the number $N$ of datapoints goes to infinity? (very large dataset)
#solution[
  In the previous question we defined $p_hat(bold(s)) (microstate(s))$, for some microstate $microstate(s)$, as the proportion of observations in the dataset $hat(bold(s))$ which were equal to $microstate(s)$. We can equivalently define $p_hat(bold(s)) (microstate(s))$ as an average over the function $bb(1)[microstate(s)^((k)) = microstate(s)]$, which is 1 if this condition holds, and 0 otherwise:
  $
  p_hat(bold(s)) (microstate(s)) = 1/N sum_(k=1)^N bb(1)[microstate(s)^((k)) = microstate(s)]
  $
  By the Law of Large Numbers, as $n -> infinity$, this average converges to its expected value:
  $
  lim_(N -> infinity) p_hat(bold(s)) (microstate(s)) = p(microstate(s))
  $
  The main result follows naturally by considering $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ as $N -> infinity$:
  $
  lim_(N -> infinity) angle.l s_i angle.r_D &= lim_(N -> infinity) lr([sum_(microstate(s) in hat(bold(s))) s_i (microstate(s)) dot p_hat(bold(s)) (microstate(s))]) 
  &= sum_(microstate(s)) s_i (microstate(s)) dot p(microstate(s)) 
  &= angle.l s_i angle.r
    $
  And analogously for the local pairwise correlation:
  $
  lim_(N -> infinity) angle.l s_i s_j angle.r_D &= lim_(N -> infinity) lr([sum_(microstate(s) in hat(bold(s))) s_i (microstate(s)) s_j (microstate(s)) dot p_hat(bold(s)) (microstate(s))])
  &= sum_(microstate(s)) s_i (microstate(s)) s_j (microstate(s)) dot p(microstate(s))
  &= angle.l s_i s_j angle.r
  $
]

== Maximum Entropy models

In many papers, the authors refer to the generalised Ising model defined in @eq:pairwise-ising-model as a _maximum entropy model_. In this section we will show that the probability distribution defined in @eq:pairwise-ising-model is indeed the (most general) probability distribution that maximises the Shannon entropy $S[p_bold(g) (bold(s))]$ given a set of constraints (which we will specify).

=== Consider a spin system with stationary probability distribution $p(bold(s))$. Can you recall the definition of the Shannon entropy $S[p(bold(s))]$? As mentioned above for the Boltzmann distribution, we will take $k_b = 1$.
  
#solution[
  With $k_b = 1$, the Shannon entropy is:
  $
  S[p(bold(s))] = -sum_(microstate(s)) p(microstate(s)) log p(microstate(s))
  $
  where the summation is over all possible microstates $bold(s)$.
]

The Ising model in @eq:pairwise-ising-model can be seen as a _Maximum Entropy Model_, constrained to reproduce the data local magnetisation and local correlation, i.e., constrained to reproduce all the data averages $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ (for all spins $s_i$ and $s_j$). We also want $p(bold(s))$ to be normalised, which introduces the additional constraint $sum_bold(s) p(bold(s)) = 1$. To summarise, we are looking for the set of $2^n$ probabilities $p(bold(s))$ such that $S[p(bold(s))]$ is maximal, and such that 
$
sum_bold(s) p(bold(s)) = 1 quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) = angle.l s_i angle.r_D quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) = angle.l s_i s_j angle.r_D
$ <eq:ising-constraints>
where $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ are constants that are computed from the data for all distinct $s_i$ and $s_j$. Note that to be more precise, we wrote $s_i (bold(s))$ (instead of just $s_i$) to specify that this is the value of $s_i$ in the state $bold(s)$ (this will help with the next questions).

=== How many constraints are there in total?
#solution[
The constraints are as follows:
//- *Maximisation:* 1 constraint,
- *Normalisation:* 1 constraint,
- *Average local magnetisation:* $n$ (one per spin)
- *Average local correlation:* $(n(n-1))/2$ (one for each pair of distinct spins) 

Thus the total number of constraints is
$
1 + n + (n(n-1))/2 = (n(n + 1))/2 + 1
$

  i.e., the length of the parameter vector $bold(g)$, with an additional normalisation constraint to impose interdependence on the elements of $bold(g)$.

  Note that _maximisation_ could be considered a constraint (though this is is a somewhat non-standard view in the context of constrained optimisation), in which case the total number of constraints increases by 1.
]

To find the shape of the distributions $p(bold(s))$ that maximises the entropy while satisfying these constraints, we introduce an auxiliary function:
$
U[p(bold(s))] = S[p(bold(s))] + lambda_0 (sum_s p(bold(s)) - 1) + sum_(i=1)^n alpha_i (sum_bold(s) p(bold(s)) s_i (bold(s)) - angle.l s_i angle.r_D) \
+ sum_("pair"(i,j))^n eta_(i j) (sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) - angle.l s_i s_j angle.r_D)
$
where we have introduced a parameter in front of each constraint we want to impose. These parameters ($lambda_0$, $alpha_i$, and $eta_(i j)$) are called Lagrange multipliers. To find $p(bold(s))$ one must maximise this auxiliary function with respect to the $2^n$ probabilities $p(bold(s))$.

=== Let us fix a choice of a state $bold(s)$. The probability $p_bold(s) = p(bold(s))$ is a parameter of $U[bold(p)]$ where $bold(p)$ is the vector of the $2^n$ probabilities. Can you show that:

  $
  (diff U[bold(p)])/(diff p_bold(s)) = -log(p_bold(s)) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)) 
  $

#solution[
  For clarity, we treat the terms in the derivative one at a time. Observe that since we are taking the partial derivative of $U$ with respect to a single element $p_bold(s)$ in the vector $bold(p)$, in each of the sums over $bold(s)'$ in $U$, all terms will be annhialated by the derivative, with the exception of $bold(s)' = bold(s)$.

  First, examining the Shannon entropy:
  $
  diff/(diff p_bold(s)) S[p(bold(s))] &= diff/(diff p_bold(s)) lr([-sum_bold(s') p(bold(s')) log p(bold(s'))]) \
  &= -log(p_bold(s)) - p_bold(s) dot 1/p_bold(s) \
  &= -log(p_bold(s)) - 1
  $

  Next, the normalisation constraint:
  $
  diff/(diff p_bold(s)) lambda_0 lr((sum_bold(s') p(bold(s')) - 1)) &= diff/(diff p_bold(s)) lambda_0 p_bold(s) &= lambda_0
  $

  The local average magnetisation constraint:
  $
  diff/(diff p_bold(s)) sum_(i=1)^n alpha_i lr((sum_bold(s') p(bold(s')) s_i (bold(s')) - angle.l s_i angle.r_D)) &= sum_(i=1)^n alpha_i lr((diff/(diff p_bold(s)) p_bold(s) s_i (bold(s)))) \
  &= sum_(i=1)^n alpha_i s_i (bold(s))
  $

  And lastly the average local correlation constraint:
  $
  diff/(diff p_bold(s)) sum_("pair"(i,j)) eta_(i j) lr((sum_bold(s') p(bold(s')) s_i (bold(s')) s_j bold(s') - angle.l s_i s_j angle.r_D)) &= sum_("pair"(i,j)) eta_(i j) lr((diff/(diff p_bold(s)) p_bold(s) s_i (bold(s)) s_j (bold(s)))) \
  &= sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))
  $

  Taking these results together, we arrive at the desired result:
  $
  (diff U[bold(p)])/(diff p_bold(s)) = -log(p_bold(s)) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))
  $
]

=== Can you show that the most general expression of $p_bold(s)$ with maximal entropy that satisfying the constraints in @eq:ising-constraints is @eq:pairwise-ising-model? Give the relation between $lambda_0$ and the partition function $Z$. How are the parameters $alpha_i$ or $eta_(i j)$ related to the parameters $h_i$ and $J_(i j)$?

#solution[
  The constrained optimisation problem described by $U$ is optimised iff, for all microstates $bold(s)$ and spins $i, j$, $ (diff U)/(diff p_s) = (diff U)/(diff lambda_0) = (diff U)/(diff alpha_i) = (diff U)/(diff eta_(i j)) = 0 $
  Note that the partial derivatives with respect to the constraints are zero exactly when those constraints are satisfied. The partial derivative with respect to $p_s$ is zero for a particular $bold(s)$ when
  $
  0 &= -log(p_bold(s)) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)) \
  p_bold(s) &= exp(-(1 - lambda_0) + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))) \
  &= 1/exp(1 - lambda_0) exp(sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)))
  $ <eq:maximum-entropy-model-in-ising-form>

  Equating @eq:pairwise-ising-model with @eq:maximum-entropy-model-in-ising-form, we obtain the following relationships:
  $
  Z &= exp(1 - lambda_0) wide&& \
  h_i &= alpha_i wide &&"for each spin" i \
  J_(i j) &= eta_(i j) wide &&"for each pair of spins" i,j
  $
]

== Statistical inference: model with no couplings

Consider the model with no couplings (all the $J_(i j) = 0$):
$
p_bold(g) (bold(s)) = 1/(Z(bold(g))) exp(sum_(i=1)^n h_i s_i)
$
The vector $bold(g) = (h_1, ..., h_n)$ now only contains $n$ local field parameters.

=== Can you show that in that case the model is assuming the variables are independent from each other, i.e., that we can write the joint probability distribution as a product of a probability distribution over each variable $p_bold(g) (bold(s)) = product_(i=1)^n p_(bold(h_i)) (s_i)$? What is the probability distribution $p_(bold(h_i)) (s_i)$ for the spin variable $s_i$?

#solution[
  We first solve for the partition function in the described model. Using the normalisation condition, we have:
  $
  Z(bold(g)) &= sum_(microstate(s)) exp(sum_(i=1)^n h_i s_i) \
  &= sum_(s_1 in plus.minus 1) dots.c sum_(s_(n-1) in plus.minus 1) exp(sum_(i=1)^(n-1) h_i s_i) lr((exp(h_n) + exp(-h_n))) \
  &= product_(i=1)^n exp(h_i) + exp(-h_i) \
  &= product_(i=1)^n 2cosh(h_i)
  $

  We can then rewrite $p_bold(g) (bold(s))$ as:
  $
  p_bold(g) (bold(s)) &= 1/(Z(bold(g))) exp(sum_(i=1)^n h_i s_i) \
  &= 1/(Z(bold(g))) product_(i=1)^n exp(h_i s_i) \
  &= product_(i=1)^n exp(h_i s_i)/(2cosh(h_i)) 
  $

  Where we recognise the term inside the product as the probability distribution for a system consisting of a single spin. It follows that $1/(2cosh(h_i))$ normalises the distribution for a single spin, and thus we find $p_bold(g) (bold(s)) = inline(product_(i=1)^n) p_bold(h_i) (s_i)$, where:

  $
  p_bold(h_i) (s_i) = exp(h_i s_i)/(2cosh(h_i))
  $
]

=== Take one of the spin variables $s_i$. We recall that $angle.l s_i angle.r_D$ is the average value of $s_i$ in the data (given a dataset, this quantity is a constant), and that $angle.l s_i angle.r = sum_bold(s) p(bold(s)) s_i$ is the model average of $s_i$. Can you show that the value of the parameter $h_i$ that satisfies the constraint $angle.l s_i angle.r = angle.l s_i angle.r_D$ is:
  $
  h_i = tanh^(-1) (angle.l s_i angle.r_D),
  $
  where $tanh^(-1) (x)$ denotes the inverse of the hyperbolic tangent? In particular, in that case the probability distribution over $s_i$ in the model is exactly equal to the empirical distribution of $s_i$.

#solution[
  Let $i$ be a spin in a model with $n$ spins, and $angle.l s_i angle.r_D$ the average value of $s_i$ in a dataset $D$. Suppose that $angle.l s_i angle.r = angle.l s_i angle.r_D$ in the model, then from the definition of $angle.l s_i angle.r$
  $
  angle.l s_i angle.r = sum_(microstate(s)) s_i (microstate(s)) dot p_bold(g) (microstate(s)) = angle.l s_i angle.r_D
  $ <eq:2-4-2-definition>
  From the previous exercise, we have that $p_bold(g) (microstate(s)) = inline(product_(j=1)^n) p_bold(h_j) (s_j)$. Taking $bold(s)_(-i) = (s_1, ..., s_(i-1), s_(i+1), ..., s_n)$ to denote the vector of spin variables excluding $s_i$, we use this result to rewrite @eq:2-4-2-definition in terms of $p_bold(h_i)$:
  $
  angle.l s_i angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) dot p_bold(h_i) (s_i (microstate(s))) dot product_(j != i) p_bold(h_j) (s_j (microstate(s))) \
  &= sum_(s in plus.minus 1) s p_bold(h_i) (s) dot overbrace(sum_(bold(s)_(-i)) product_(j != i) p_bold(h_j) (s_j (bold(s)_(-i))), = 1) \
  &= p_bold(h_i) (1) - p_bold(h_i) (-1) \
  &= (exp(h_i) - exp(-h_i))/(2cosh(h_i)) \
  &= sinh(h_i)/cosh(h_i) \
  &= tanh(h_i)
  $

  In the second line of reasoning we group the original summation terms by their $s_i$ component, conditional on observing any other combination of states in $bold(s)_(-i)$. Since each $p_bold(h_j)$ is an independent probability distribution, it follows that the sum over the joint probabilities of $bold(s)_(-i)$ is normalised, hence allowing the simplification in the third line.

  Finally, we obtain the desired result, $h_i = tanh^(-1)(angle.l s_i angle.r_D)$. 
]

=== In Eq. (6), we observe that:
  - If $angle.l s_i angle.r_D > 0$, then the inferred $h_i$ is also positive;
  - Reciprocally, $angle.l s_i angle.r_D < 0$, then the inferred $h_i$ is also negative.
  How does this connect with the tendency of the $i$'th judge to vote on average more liberal or more conservative? Is this result coherent with the general comments that we did in Question Q1.4.?

#solution[
  This result is coherent with our comments in Q1.4, in which we interpreted the sign of $h_i$ as reflecting $i$'s political leaning ($-1$ for liberal, $+1$ for conservative). As $angle.l s_i angle.r_D$ is the average value of $s_i$ in a dataset $D$, a positive value occurs when $i$ votes conservative in more than 50% of cases. Likewise, a negative value occurs when $i$ votes liberal in more than 50% of cases. 

  As in Q1.4, if the $i$'th judge has an equal number of liberal and conservative votes in $D$, then $angle.l s_i angle.r_D = h_i = 0$. 

  While in Q1.4 our discussion was concerned with a model which included interactions, the comments are still relevant here, as we interpreted the sign of $h_i$ to reflect political leaning in absence of interactions with other judges. In this model we have made this assumption explicit by taking $J_(i j) = 0$.
]

== Statistical inference: maximising the log-likelihood function

*Introducing the likelihood function.* Looking more closely at @eq:pairwise-ising-model, one can see that it does not just define a single probability distribution, but many of them: there is one probability distribution for each value of the set of parameters $bold(g)$. More precisely, the distribution in @eq:pairwise-ising-model changes continuously as one continuously varies the parameters in $bold(g)$. We say that @eq:pairwise-ising-model defines a _parametric family of probability distributions_. The inference procedure consists in finding the value of the parameters $bold(g)$ that maximises the probability that the model $p_bold(g) (bold(s))$ produces the data.

To do so, we introduce the _log-likelihood function_:
$
cal(L)(bold(g)) = log P_bold(g) (microstate(s))
$ <eq:2-5-log-likelihood>
where $P_bold(g) (microstate(s))$ is the probability that the model $p_bold(g) (bold(s))$ produces the dataset $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$. Note that $cal(L)(bold(g))$ is a function of the parameters $bold(g)$. The inference procedure therefore consists in finding the value $bold(g)^star$ of the parameters that maximises $cal(L)(bold(g))$. For the moment we will assume that there exists only a unique such value of $bold(g)$.

=== We assume that, in the dataset $hat(bold(s))$, all datapoints are independently sampled from an underlying distribution $p_bold(g) (bold(s))$. Can you show that the log-likelihood function can be re-written as:
  $ cal(L)(bold(g)) = N sum_(bold(s)) p_D (bold(s)) log p_bold(g) (bold(s)) $ <eq:log-likelihood-computable>
  where $p_D (bold(s))$ is the empirical distribution over the states? The empirical distribution is given by $p_D (bold(s)) = K(bold(s))/N$ where $K(bold(s))$ is the number of times that the datapoint $bold(s)$ occurs in the dataset.

#solution[
  From the independent sampling assumption, we have that $P_bold(g) (hat(bold(s))) = inline(product_(k=1)^N) p_bold(g) (hat(bold(s))^((k)))$. Substituting this in @eq:2-5-log-likelihood, we have:
  $
  cal(L)(bold(g)) &= log(product_(k=1)^N p_bold(g) (hat(bold(s))^((k)))) \
  &= sum_(k=1)^N log p_bold(g) (hat(bold(s))^((k))) \
  &= sum_(bold(s)) K(bold(s)) log p_bold(g) (bold(s)) \
  &= sum_(bold(s)) N dot K(bold(s))/N log p_bold(g) (bold(s))\
  &= N sum_(bold(s)) p_D (bold(s)) log p_bold(g) (bold(s))
  $
  as desired, where the summation is over uniques microstates in the dataset (or equivalently over all possible microstates).
]

  *Ising model.* We now take the model distribution $p_bold(g) (bold(s))$ to be given by the Ising model in @eq:pairwise-ising-model.
=== Taking the first derivative of $cal(L)(bold(g))$ with respect to a parameter $h_i$, can you show that at the maximum of $cal(L)(bold(g))$ we have that $angle.l s_i angle.r = angle.l s_i angle.r_D$? Similarly, taking the first derivative of $cal(L)(bold(g))$ with respect to a parameter $J_(i j)$, can you show that at the maximum of $cal(L)(bold(g))$ we have that $angle.l s_i s_j angle.r = angle.l s_i s_j angle.r_D$?

#solution[
  Let $p_bold(g) (bold(s))$ be as defined in @eq:pairwise-ising-model, and let $i$ correspond to the spin $s_i$. Then $cal(L)(bold(g))$ is:
  $
  cal(L)(bold(g)) &= N sum_(bold(s)) p_D (bold(s)) log lr([1/(Z(bold(g))) exp(sum_(j=1)^n h_j s_j + sum_("pair"(j,k)) J_(j k) s_j s_k)]) \
  &= N sum_(bold(s)) p_D (bold(s)) lr([sum_(j=1)^n h_j s_j + sum_("pair"(j,k)) J_(j k) s_j s_k - log Z(bold(g))])
  $

  The partial derivative of $cal(L)(bold(g))$ with respect to $h_i$ is then:
  $
  (diff cal(L))/(diff h_i) &= N sum_(bold(s)) p_D (bold(s)) lr([diff/(diff h_i) sum_(j=1)^n h_j s_j + diff/(diff h_i) sum_("pair"(j,k)) J_(j k) s_j s_k - diff/(diff h_i) log Z(bold(g))]) \
  &= N sum_(bold(s)) p_D (bold(s)) dot (s_i - diff/(diff h_i) log Z(bold(g)))
  $ <eq:2-5-deriv-log-likelihood-unsimplified>

  Where the partial derivative of the log partition function with respect to $h_i$ is:
  $
  (diff log Z(bold(g)))/(diff h_i) &= 1/(Z(bold(g))) dot diff/(diff h_i) sum_(bold(s)') exp(sum_(j = 1)^n h_j s_j + sum_("pair"(j,k)) J_(j k) s_j s_k) \
  &= 1/(Z(bold(g))) dot sum_(bold(s)') diff/(diff h_i) exp(sum_(j = 1)^n h_j s_j + sum_("pair"(j,k)) J_(j k) s_j s_k) \
  &= 1/(Z(bold(g))) dot sum_(bold(s)') exp(sum_(j = 1)^n h_j s_j + sum_("pair"(j,k)) J_(j k) s_j s_k) dot s_i \
  &= sum_(bold(s)') p_bold(g) (bold(s')) s_i \
  &= angle.l s_i angle.r 
  $ <eq:2-5-deriv-log-partition>

  $cal(L)(bold(g))$ attains its maximum value with respect to $h_i$ when $(diff cal(L))/(diff h_i) = 0$. Substituting @eq:2-5-deriv-log-partition into @eq:2-5-deriv-log-likelihood-unsimplified and solving for the maximum, we find:
  $
  0 &= N sum_(bold(s)) p_D (bold(s)) dot (s_i - angle.l s_i angle.r) \
  ==> wide N sum_(bold(s)) p_D (bold(s)) angle.l s_i angle.r &= N sum_(bold(s)) p_D (bold(s)) s_i \
  angle.l s_i angle.r dot sum_(bold(s)) p_D (bold(s)) &= angle.l s_i angle.r_D \
  angle.l s_i angle.r &= angle.l s_i angle.r_D
  $
  Where the left-hand side summation cancels in the final step since the occurrence proportions of states in $D$ must sum to 1.

  We proceed analogously to show that $angle.l s_i s_j angle.r = angle.l s_i s_j angle.r_D$ when $cal(L)(bold(g))$ is maximised with respect to $J_(i j)$:

  $
  (diff cal(L))/(diff J_(i j)) &= N sum_(bold(s)) p_D (bold(s)) lr([diff/(diff J_(i j)) sum_(k=1)^n h_k s_k + diff/(diff J_(i j)) sum_("pair"(k,ell)) J_(k ell) s_k s_ell - diff/(diff J_(i j)) log Z(bold(g))]) \
  &= N sum_(bold(s)) p_D (bold(s)) dot (s_i - diff/(diff h_i) log Z(bold(g))) \
  &= N sum_(bold(s)) p_D (bold(s)) dot (s_i - sum_(bold(s)') p_bold(g) (bold(s)') s_i s_j) \
  &= N sum_(bold(s)) p_D (bold(s)) dot (s_i - angle.l s_i s_j angle.r) \
  $ <eq:2-5-deriv-local-correlation>

  As before, evaluating @eq:2-5-deriv-local-correlation at $(diff cal(L))/(diff J_(i j)) = 0$ gives:
  $
  N sum_(bold(s)) p_D (bold(s)) angle.l s_i s_j angle.r &= N sum_(bold(s)) p_D (bold(s)) s_i s_j \
  angle.l s_i s_j angle.r &= angle.l s_i s_j angle.r_D
  $
]

= Application to the analysis of the US supreme Court
A system with $n$ spin variables can be in $2^n$ different states. However, most of the time, the number of different states observed in a dataset is very small compared to $2^n$.

#let supreme_court_results = json("../results/supreme_court.json")
#let n_spins = supreme_court_results.n_spins
#let n_data = supreme_court_results.n_data
#let n_unique_data = supreme_court_results.n_unique_obs
#let possible_vote_combos = calc.pow(2, n_spins)

=== For the US Supreme court dataset: What is the number $n$ of spin variables, and the total number $2^n$ of states that can be observed for that system? What is the total number $N$ of datapoints in the dataset? What is the number $N_"max"$ of different states that are observed?

#solution[
  The number of spins corresponds to the number of judges ($n=#n_spins$), for which there are $2^#n_spins = #possible_vote_combos$ possible combinations of votes. The dataset contains $#n_data > #possible_vote_combos$ rows; however, only $N_"max" = #n_unique_data$ of these are unique. 
]

=== *(Bonus question)* Numerical solution: For the dataset provided, find numerically the value of the parameters $bold(g)$ of the fully-connected pairwise model @eq:pairwise-ising-model that maximises the log-likelihood function $log cal(L)(bold(g))$. What are the main computational limitations of your algorithm? How can you improve it?

#solution[
  We solve for the parameters $bold(g)$ by minimising the negative log-likelihood function @eq:log-likelihood-computable using the SciPy BFGS solver. Note that this is equivalent to _maximising_ the log-likelihood function. The negative log-likelihood across optimisation iterations is shown in @fig:model-fitting-neg-log-likelihood.

  #figure(image("../results/figures/neg_log_likelihood_model_fitting.svg")) <fig:model-fitting-neg-log-likelihood>

  While we obtain similar fit parameters to those provided on Canvas, the optimisation process is computationally expensive -- even for only nine judges -- in particular because the partition function must be re-evaluated at each step. Computing $Z(g)$ requires a summation over $2^n$ distinct microstates, with $O(n^2)$ operations for each term.

  One approach to alleviating this computational limitation is to approximate the distribution $p_bold(g) (bold(s))$ using Metropolis sampling, such that we don't need to calculate the partition function. This dataset contains a small number of states which occur with relatively high probability -- the two most common microstates absorb 45% of the probability mass. As such, there is reason to expect a reasonably accurate estimate for $p_bold(g)(bold(s))$ from a small number of samples (compared to the number of microstates).
]

=== We would like to reproduce *Figure 13* of the attached paper. Can you plot the $angle.l s_i angle.r_D$ as a function of $i$? Can you re-order the label $i$ so that the values of $angle.l s_i angle.r_D$ are ordered from the smallest (negative value) to the largest (positive value), as in Fig. 13.A top? Keeping the new ordering, can you plot a heatmap of the matrix of $angle.l s_i s_j angle.r_D$ (see Fig. 13.A bottom)? What can we say about the judges with negative $angle.l s_i angle.r_D$? With positive $angle.l s_i angle.r_D$? Can you comment on these plots?

#solution[
  We calculate $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ using the methods derived in @eq:local-magnetisation-using-pdata and @eq:local-correlation-using-pdata. The results, in increasing order of $angle.l s_i angle.r_D$, are displayed in @fig:avg-vote-per-judge and @fig:avg-pairwise-correlation respectively. 

  In @fig:avg-vote-per-judge we observe a 5-4 split between judges who vote mostly conservative, and those who mostly vote liberal. We also note that $abs(angle.l s_i angle.r_D) < 0.5$ for all judges, indicating that even the most conservative/liberal judges vote against their ideological leaning in at least 25% of cases. 


  #figure(
    image("../results/figures/empirical_avg_local_mag_sorted.svg"),
    caption: [Average vote per judge in the supplied dataset, with x-axis ordered according to increasing tendency to vote conservative. Votes in the dataset are assigned $-1$ if considered *liberal*, and $+1$ if considered *conservative*.]
  ) <fig:avg-vote-per-judge>

  @fig:avg-pairwise-correlation also exhibits several features of interest. The liberal and conservative blocks are clearly visible; however, no two judges are completely aligned. Within the liberal block, the first judge (who votes liberal most often) is distinct from the other three. A similar separation appears to exist in the conservative block as well, between the first three and final two judges -- this is less evident in the context of the conservative block alone, yet we see the final two judges show extremely low correlation with the most liberal judge. 

  A particularly counter-intuitive feature of @fig:avg-pairwise-correlation is that all observed values are positive. We observed in @fig:avg-vote-per-judge that some judges predominantly vote liberally, and others conservatively. If we suppose that the judges' votes are independent of one another, we would expect to see a negative correlation between liberal-conservative pairs of judges. Since we observe _no negative correlations_, this implies that the votes are not independent, and that interaction between judges may play a significant role in any particular judge's vote -- irrespective of their political leaning. 

  #figure(
    image("../results/figures/empirical_avg_local_corr_sorted_by_mag.svg"),
    caption: [Pairwise correlations between judges' votes. Axes ordered by increasing average tendency to vote conservative.]
  ) <fig:avg-pairwise-correlation>

  Critically, we note that the pairwise correlations in @fig:avg-pairwise-correlation are not sufficient for extended discussion on the interactions between judges, as the values include the heterogeneous influence of individual judges' voting patterns.
]

=== Keeping the new ordering of the labels $i$, can you plot a heatmap of the fitted vector of $h_i$'s and a heatmap of the fitted matrix of $J_(i j)$'s, as in Fig. 13.B? You can use the fitted values of $h_i$ and $J_(i j)$ provided in Canvas. Can you comment on these plots?
  
  Note that the values of $h_i$ and $J_(i j)$ that are provided in Canvas are following the same order of the variables $s_i$ than in the USSC datafile (i.e., the original ordering of the judges).

#solution[
  The fit $h_i$ parameters are displayed in @fig:fit-h-scatterplot, and the interaction terms $J_(i j)$ in @fig:fit-j-heatmap. In both cases, the judges are sorted in their axes in accordance with increasing order of $angle.l s_i angle.r_D$, the tendency for a judge to vote conservatively. We observe substantial differences between each figure and its counterpart in the previous question.

  The $h_i$ values displayed in @fig:fit-h-scatterplot reflect the political leaning of each judge, in absence of interactions. We observe several surprising inversions with respect to $angle.l s_i angle.r_D$. For instance, @fig:avg-vote-per-judge displayed a tendency for the third judge to vote (on average) more liberally, yet $h_3 > 0$ indicates a conservative leaning. On the other hand, $h_5$ and $h_7$ are both close to 0, yet these judges' votes displayed non-neglibible conservative tendencies. 

  To understand this counter-intuitive result, we must consider it in the context of the model interaction terms. Indeed, in *Q1.4* we interpreted $h_i$ as the political leaning of the $i$'th judge, _in absence of interactions_.

  //is positive. between judges compared to theJudges 3, 5, and 7 are of particular interest. In @fig:avg-vote-per-judge, we observed that the third judge voted liberal more frequently than conservative, yet $h_3 > 0$ indicates a conservative bias after accounting for pairwise interactions. On the other hand, $h_5$ and $h_7$ are much lower than expected (though still greater than 0).
  #figure(
    image("../results/figures/fit_parameter_h_scatterplot.svg"),
    caption: [Fit $h_i$ parameter values per judge in the supplied dataset, with x-axis sorted according to increasing tendency to vote conservative. Votes in the dataset are assigned $-1$ if considered *liberal* and $+1$ if considered *conservative*.]
  ) <fig:fit-h-scatterplot>

  The fit interaction terms are displayed in @fig:fit-j-heatmap. The values $J_(i j)$ indicate the tendency for a pair of judges to vote in-line or against one another, regardless of their political leaning. A value close to zero implies that any similarity/dissimilarity in the votes of a pair of judges is well-explained by their political leanings (or interactions with other judges). A positive value suggests that two judges may vote the same way often, even when this contradicts ones' political leaning. Note that the diagonal terms are 0 on account of the model disregarding self-interactions (i.e., such parameters do not exist in the model). 



  #figure(
    image("../results/figures/fit_parameter_j_heatmap.svg"),
    caption: [Interaction terms $J_(i j)$ between judges in a pairwise Ising model fit to the supplied dataset. Axes are ordered according to increasing tendency for a judge to vote conservatively.]
  ) <fig:fit-j-heatmap>

  We observe qualitatively different behaviour with respect to the local correlations in @fig:avg-pairwise-correlation. Firstly, the heatmap now includes negative terms, though these typically have lower magnitude that the positive interactions. The liberal and conservative "blocks" identified earlier each display positive pairwise interaction values, with negative interactions restricted to cross-block pairs of judges. Interestingly, we also observe some positive interactions between the blocks, which vary between judges. For instance, the $0$'th judge has a positive interaction value with judges 5 and 8, while the $1$'st judge displays a negative interaction with the 8'th judge, but positive interactions with the 6'th and 7'th. Overall we observe considerably more heterogenerous behaviour in @fig:fit-j-heatmap than in @fig:avg-pairwise-correlation.

  Returning to our examination of judges 3, 5, and 7, @fig:fit-j-heatmap offers some insight into their unexpected inferred political leanings. Consider the 3rd judge. In the data, they display a tendency to vote liberally, yet $h_3 > 0$ indicates that this tendency does not arise due to their political leaning. For the model to accurately reflect the data (we will see in later discussion that it does), this behaviour must then arise from the 3rd judge's interactions with other judges. 

  Indeed, we see from @fig:fit-j-heatmap that they exhibit positive interactions with judges 0, 1 and 2, and negative interaction values with judges 7 and 8 (who we earlier predicted to be the most conservative). Intuitively, this implies that the 3rd judge tends to vote similarly to liberal judges, but differently to the more conservative judges. Thus their liberal voting tendencies are explained not by their political leaning, but by their tendency to vote "like a liberal judge". 

  The 5'th judge actually exhibits positive interaction values with all judges, yet their strongest interaction is with the 6'th judge, who in turn consistently votes opposite to the 0'th judge (the most liberal), explaining the tendency for the 5'th judge to vote conservative. Likewise, the 7'th judge -- who we had predicted as one of the more conservative -- displays only a small conservative leaning. This is explained by their strong positive interactions with the 8'th judge, and tendency to vote oppositely to the 0'th judge. 
]

=== *Scatter plot 1, "Cross-validation":* For all the states observed in the data, can you plot the empirical probability of the state, $p_D (bold(s))$ against the model probability $p_bold(g) (bold(s))$ with the fitted parameters ($h_i$ and $J_(i j)$)? What can you say about this plot?

#solution[
  @fig:cross-validation compares the probability of observed microstates under the fit Ising model with that observed in the dataset. We see general agreement between the two. In particular, the most probable states in the dataset have similar probability under the fit model. 

  The discrepancies indicate potential higher-order interactions which are not captured by the pairwise Ising model. Since this model is constrained to reproduce $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ for each $s_i, s_j$, these discrepancies cannot be explained by the political leaning of particular judges or pairwise interactions alone. However, such discrepancies appear relatively minor, suggesting that the pairwise model is a reasonably good explanation for the data.

  #figure(
    image("../results/figures/empirical_vs_model_probability.svg"),
    caption: [Cross-validation check: comparison between microstate probability under the fit pairwise Ising model, and the observed data. The $y=x$ line is indicated by a grey dashed line.]
  ) <fig:cross-validation>
]

=== *Scatter plot 2: Checking the fit:* For each spin $s_i$, can you plot the value of $angle.l s_i angle.r_D$ in the data against the value $angle.l s_i angle.r$ in the fitted model? For each pair of spins $s_i, s_j$, can you plot the value of $angle.l s_i s_j angle.r_D$ in the data against the value of $angle.l s_i s_j angle.r$ in the fitted model? What can you say about these plots?

#solution[
  @fig:comparison-avg-local-mag-data-model and @fig:comparison-avg-local-pairwise-corr-data-model compare the model average magnetisation and pairwise correlation with those values calculated from the observed data. We observe strong agreement between the two. Following our derivation in *Q5.2*, this implies that the model parameters are well-fit to the data, and in particular that the log-likelihood is maximal.
  #figure(
    image("../results/figures/avg_local_magnetisation_data_vs_model.svg"),
    caption: [Comparison of average local magnetisation under the fit model and as calculated from the observed data. The $y=x$ line is indicated in grey. We observe that the model accurately reproduces the data value.]
  ) <fig:comparison-avg-local-mag-data-model>

  #figure(
    image("../results/figures/avg_local_correlation_data_vs_model.svg"),
    caption: [Comparison of average local pairwise correlation under the fit model and as calculated from the observed data. The $y=x$ line is indicated in grey. We observe that the model accurately reproduces the data value.]
  ) <fig:comparison-avg-local-pairwise-corr-data-model>
]

  One of the questions addressed in the paper is: can the pairwise model reproduce higher-order patterns of the data better than a model of independent judges? To answer this question, the authors introduce the probability $P(k)$ that there are $k$-conservative votes as answer to a case (i.e., the probability that a datapoint contains $k$ conservative votes).
=== Consider a model with no coupling, in which judges vote independently from each other. Each judge $s_i$ has a probability $p_i (+1)$ to vote conservative. In that case, what is the probability $P_I (k)$ that $k$ judges vote conservative? Note: we added the label "I" to $P(k)$ to specify that it is the probability distribution obtained for an independent model.
  
  In the dataset, how many judges have votes that are more conservative on average? At which value of $k$ do you then expect the maximum of $P_I (k)$ to be? Can you obtain the values of $p_i (+1)$ from the data and plot the values of $P_I (k)$ as a function of $k$ for the independent model?

#solution[
  From @fig:avg-vote-per-judge there are 5 judges whose votes are, on average, more conservative. Since $P_I (k)$ assumes that judges vote independently, we then expect a maximum at $k=5$. This makes intuitive sense, since $k=6$ requires some judges who vote liberally to vote conservatively, and $k=4$ would require an (on-average) conservative judge to vote liberally. These are both less likely than the case where each judge votes in-line with their typical behaviour. Note that the peak at $k=5$ does include mass from scenarios where typically-liberal judges vote conservative, and vice versa.

  It is straightforward to obtain $p_i (+1)$ from the data as the proportion of rows where each judge votes conservative. It is less straightforward to compute $P_I (k)$ given each $p_i (+1)$, as $k$ conservative votes could come from any one of $binom(9,k)$ combinations of judges.

  One approach to simplifying the calculation of $P_I (k)$ is to redefine it as a recurrence relation which can be computed straightforwardly using recursion or dynamic programming. We define the probability $q(x,i)$ as the probability of observing $x$ conservative votes among the final $i$ judges (under an arbitrary ordering of judges). Formally, we define $q$ by:

  $
  q(x, i) = cases(
    0 &"if" x > i\
    product_(j=(n-i)+1)^n p_j (-1) &"if" x = 0 \
    product_(j=(n-i)+1)^n p_j (+1) &"if" x = i \
    p_((n-i) + 1) (+1) dot q(x-1, i-1) + p_((n-i) + 1) (-1) dot q(x, i-1) quad&"if" x < i
  )
  $ <eq:probability-of-x-cons-votes-in-i-judges>

  In the first case ($x > i$), we require more conservative votes than remaining judges, which naturally has probability zero. In the second ($x=0$), any remaining votes must be liberal, as we have achieved the required number of conservative votes. Likewise, in the third case there are exactly as many judges remaining as required conservative votes, so all must vote conservative. 

  Finally, in the last case, we consider two scenarios which are recursively defined: where the current judge votes conservative, requiring one fewer conservative votes among the remaining judges, or where the current judge votes liberal, and we still require $x$ conservative votes.

  To avoid possible floating-point innacuracies when multiplying small probabilities, the products are implemented as the exponentiated sum of log-probabilities.

  The distribution $P_I (k) = q(k, n)$ is shown below in @fig:conservative-vote-count-distribution-independent.

  #figure(
    image("../results/figures/conservative_vote_count_distribution.svg"),
    caption: [Distribution of the number of conservative votes $k$, calculated with respect to the (assumed independent) probabilities of each individual judge voting conservative.]
  ) <fig:conservative-vote-count-distribution-independent>
]

=== Let us call $P_D (k)$ the probability distribution $P(k)$ obtained directly from the data. How can you compute the values of $P_D (k)$ from the data for $k=1$? for $k=2$? for any $k$? Can you plot $P_D (k)$ as a function of $k$ and compare the curve to the one obtained previously for the independent model? Where is the maximum of $P_D (k)$? Is the independent model a good model for the US Supreme Court data?

#solution[
  Observe that in a model with $n$ judges, if $k <= n$ judges vote conservative, then with the ${-1, +1}$ valuation on votes, the sum of the votes gives a net vote of 
  $
  1 dot k + (-1) dot (n - k) = 2k - n
  $

  It is then straightforward to calculate $P_D (k)$ for any $k$ as the proportion of rows in $D$ whose values sum to $2k - n$. The resulting distribution is shown with comparison to $P_I (k)$ in @fig:conservative-vote-count-distribution-independent-and-data.

  #figure(
    image("../results/figures/conservative_vote_count_distribution_2.svg"),
    caption: [Comparison between distribution of microstates in the data (orange) with the distribution predicted under the assumption that judges' votes are independent of one another (blue).]
  ) <fig:conservative-vote-count-distribution-independent-and-data>

  We observe considerable differences between the two distributions. While $P_D (k)$ still has a local maximum at $k=5$, this state (five conservative votes) is significantly less likely than predicted by the independent model. Additionally, $P_D (k)$ exhibits two other local maxima at $k=0$ and $k=9$, with the latter also being the global maximum. 

  A plausible (but incorrect) interpretation of these maxima is that they correspond to issues which are "too liberal" or "too conservative" for all judges. However, given the substantial probability mass at these peaks, if this were the case we would not observe such low probability for these cases in $P_I (k)$. 

  More likely, it is the case that these peaks correspond to issues where judges' votes cannot be explained entirely by their political leaning, i.e., that interactions between judges has resulted in agreement, sufficient to cause liberal judges to vote conservative, or vice versa. For this reason (as well as the evident qualitative differences) we conclude that $P_I (k)$ is not a good model for this data.
]

=== Let us call $P_P (k)$ the probability distribution $P(k)$ obtained from the fitted Ising model with pairwise couplings (i.e., with the model $p_bold(g) (bold(s))$ in @eq:pairwise-ising-model with the fitted parameters $h_i$ and $J_(i j)$ provided in Canvas). How can you compute $P_P (k)$: can you write the analytical expression of $P_P (k)$ as a function of $p_bold(g) (bold(s))$? Can you plot $P_P (k)$ as a function of $k$ and compare the curve to $P_D (k)$ and $P_I (k)$ previously obtained? (see Figure 16.A of the paper). Which conclusions can you draw?

#solution[
  @fig:cons-vote-count-all-compare compares the distributions of observed microstates in the data (orange) with the distributions obtained under the independent-judges assumption (blue) and the fit pairwise Ising model (green). We immediately observe that the pairwise Ising model fits the expected distibution much better than the independent-voting model. Unlike the independent voting model, the Ising model manages to accurately capture the observed peaks at $k=0$ and $k=9$, as well as the local (but smaller) maxima at $k=5$, and the troughs between these regions.


  #figure(
    image("../results/figures/conservative_vote_count_distribution_3.svg"),
    caption: [Comparison between the distribution of microstates observed in the data (orange) with that of the fit pairwise Ising model (green) and a model which assumes independent votes between judges (blue).]
  ) <fig:cons-vote-count-all-compare>

  The three peaks reflect three "stable" states in this system. In a typical case we expect to see either broad agreement (consensus) at one of the two peaks, or a split vote along the political boundary. There is a substantial decrease in probability mass at $k=1$ and $k=8$, indicating that highly-skewed, non-consensus states are uncommon -- if almost all judges are voting in a particular direction, the final judge tends to vote this way as well. Of the three stable states, roughly 45% of the probability mass lies at the extremes, indicating that such agreement is quite typical, even if it contradicts individual preferences.

  It bears to comment further on the high probability found at the extremes ($k=0$ and $k=9$). Under the independent-voting model these states have almost no probability of occurring. To see them attain such high likelihood under the Ising model is evidence of the high impact of pairwise interactions in this system.

  Finally, we again find that the pairwise Ising model, while close, does not precisely match the distribution observed in the data, particularly at $k=1$ and $k=8$. In both cases the pairwise Ising model overestimates the likelihood of observing such states. While this may simply be a facet of the dataset size, or ambiguity in the labelling, it may also indicate the presence of higher-order interactions which are not being captured by our model.
]

