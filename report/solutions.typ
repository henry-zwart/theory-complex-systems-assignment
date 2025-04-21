#import "@preview/cetz:0.3.4"
#import "@preview/fletcher:0.5.7" as fletcher: diagram, node, edge, shapes
#import "@local/uva-game-theory:0.1.0": *

#set math.equation(numbering: "(1)", number-align: bottom)
#set block(width: 100%)

#let microstate(x) = $hat(bold(#x))$
#let ensembleavg(varname, microstatevar: $s$) = $sum_(microstate(microstatevar)) varname(microstate(microstatevar)) dot P(microstate(microstatevar))$

= Modelling the activity of a single neuron

#let neuron_results = json("../results/neuron_modelling.json")

+ Can you plot the distribution $P(tau)$ of the time intervals $tau$ between successive spikes? Check that there is indeed a rafactory period, i.e., a time interval $tau_0$ after each spike, during which the neuron doesn't spike again. What is the duration $tau_0$ for this time interval?
  #solution[
  #figure(
    image("../results/figures/neuron_interspiking_distribution.svg"),
    caption: [Distribution over the duration between activity spikes for a neuron. A refractory period of $tau_0 approx 1.9 m s$ is identified as the minimum observed duration between spikes.]
  ) <fig:neuron-empirical-interspike-distribution>

  Refractory period duration $tau_0 = #calc.round(neuron_results.tau0, digits: 2)"ms"$ calculated as the minimum observed duration between neuron spikes.
  ]

+ Can you check that the decay of the distribution $P(tau)$ of inter-spike intervals is indeed exponential? Measure the corresponding decay rate $lambda$
  #solution[
  To examine the decay of $P(tau)$, we consider $tau > tau_0 approx #calc.round(neuron_results.tau0, digits: 3)$. Exponential decay is a good model for $P(tau - tau_0)$ if we observe a linear relationship in @eq:exponential-decay-linear-form, where the decay rate $lambda$ is given by the slope:
  #{
  set math.equation(number-align: bottom)
    [$
    P(tau - tau_0) &= P(tau_0) e^(-lambda (tau - tau_0)) \
    1/(P(tau - tau_0)) &= 1/P(tau_0) e^(lambda (tau - tau_0)) \
    log(1/(P(tau - tau_0))) &= log(1/P(tau_0)) + log(e^(lambda (tau - tau_0))) \
    -log(P(tau - tau_0)) &= -log(P(tau_0)) + lambda (tau - tau_0) \
    $ <eq:exponential-decay-linear-form>]
  }

  
  To estimate the probability distribution over $tau$, we subtract the refractory period from each inter-spike duration, and bin the data into 50 bins of uniform width. We normalise the number of observations in each bin with respect to the total number of observations ($n=#neuron_results.n$) and take this value to be the empirical probability for the value of $tau$ at the center of each bin.   
  #figure(
    image("../results/figures/neuron_exponential_fit.svg"),
    caption: [Neuron inter-spike durations (without refractory period) shows a linear relationship under a log-transform, with slope $lambda = #calc.round(neuron_results.exponential_fit.lambda, digits: 4)$ and intercept $P(tau_0) = #calc.round(-calc.ln(neuron_results.exponential_fit.p_tau_0), digits: 4)$. $R=#calc.round(neuron_results.exponential_fit.r_value, digits: 4)$.]
  ) <fig:neuron-log-linear-fit>

  @fig:neuron-log-linear-fit shows the empirical probabilities transformed using @eq:exponential-decay-linear-form. We observe a good linear fit, indicating exponential decay in $P(tau)$ with $lambda = #calc.round(neuron_results.exponential_fit.lambda, digits: 4)$ given by the slope of the linear fit.

  The data contains fewer samples for large $tau$, which explains the reduced quality of the exponential fit.
  ]

+ Can you deduce an analytical expression for the distribution of inter-spike time interval $P(tau)$ of the delayed Poisson process as a function of $lambda$ and $tau_0$? Compare your model distribution to the one obtained from the data.
  #solution[
  $
  P(tau) = cases(
    0 quad &tau < tau_0 \
    lambda exp(-lambda (tau - tau_0)) quad &tau >= tau_0
  )
  $
  ]

+ Using your model, can you generate another 1000 (spike times) datapoints?
  #solution[
  #figure(image("../results/figures/neuron_compare_empirical_theoretical.svg"))
]

+ What is the average spiking rate $f$ of the neuron in the data? How is $f$ analytically related to $tau_0$ and $lambda$ that you have previously measured?
  #solution[
  The average spiking rate, $f = #calc.round(1000 * neuron_results.mean_spike_rate, digits: 2)"hz"$ is calculated as the mean $tau^(-1)$ over the dataset. We can express $f$ analytically in terms of the expected value, as $f = 1\/E[tau]$. First, solving for $E[tau]$:
  $
  E[tau] &= integral_0^infinity tau P(tau) dif tau \
  &= integral_0^(tau_0) tau P(tau) dif tau + integral_(tau_0)^infinity tau P(tau) dif tau \
  &= integral_(tau_0)^infinity tau dot lambda e^(-lambda (tau - tau_0)) dif tau \
  &= [-tau e^(-lambda (tau - tau_0))]_(tau_0)^infinity - integral_(tau_0)^infinity -e^(-lambda (tau - tau_0)) dif tau \
  &= lr([-tau e^(-lambda (tau - tau_0))] - (1/lambda) e^(-lambda (tau - tau_0))])_(tau_0)^infinity \
  &= tau_0 + 1/lambda
  $

  The derived value for $E[tau]$ is as expected, since it represents the expected value of a regular exponential distribution with rate $lambda$, shifted by the refractory period $tau_0$. The expected spiking rate is then:
  $
  f &= 1/(E[tau]) \
  &= 1/(tau_0 + 1\/lambda) \
  &= lambda/(lambda tau + 1)
  $ <eq:neuron-expected-spiking-rate>

  Using the inferred values for $lambda$ and $tau_0$ from the prior question, @eq:neuron-expected-spiking-rate gives $f=#{calc.round(1000 * neuron_results.exponential_fit.lambda/(neuron_results.exponential_fit.lambda * neuron_results.tau0 + 1), digits: 2)}"hz"$. 
  ]

= Modelling binary data with the Ising model

== Pairwise spin model

+ How many terms are in the sum over the $"pair"(i, j)$? Can you deduce what is the number of parameters in the vector $g = (h_1,..., h_n, J_(1,2),..., J_n$? Can you re-write the sum over the $"pair"(i, j)$ as a double sum over $i$ and $j$ (without counting twice each pair)?
  #solution[
  The number of pairs in the sum over $"pair"(i, j)$ is given by the number of ways which we can choose two distinct spins from a set of $n$ total spins, such that the order of the chosen spins doesn't matter:

  $
  (n (n - 1))/2
  $

  So the total number of parameters in $bold(g)$ is:

  $
  n + (n (n-1))/2 = (n(n+1))/2
  $

  We can rewrite the sum using a double summation over $i$ and $j$ by enumerating the $k = n-1$ ways to choose the first spin, and then the $n - k$ ways to choose the second:

  $
  sum_("pair"(i,j)) J_(i,j) s_i s_j = sum_(i=1)^(n-1) sum_(j=i+1)^n J_(i,j) s_i s_j
  $
  ]

+ Can you write down explicitly the terms in the exponential of Eq. (1) for a system with $n = 3$ spins?
  #solution[
  $
  h_1 s_1 + h_2 s_2 + h_3 s_3 + J_(1,2) s_1 s_2 + J_(1,3) s_1 s_3 + J_(2,3) s_2 s_3
  $
  ]

+ In Eq. (1), we can recognize the Boltzmann distribution, in which the parameter $beta = 1\/(k_b T)$ was taken equal to 1 (more precisely, the constant $k_B$ was taken equal to 1, and the temperature parameter $T$ was absorbed in the parameters $h_i$ and $J_(i j)$). What is the energy function associated with the Boltzmann distribution in that case? What is the partition function and what is its general expression?
  #solution[
  The Boltzmann distribution has the form $P(hat(bold(s))) = 1/Z exp(-beta E(hat(bold(s))))$, where $beta = 1/(k_b T)$. Then from Eq. (1) we have
  $
  -1/(k_b T) E(hat(bold(s))) &= sum_(i=1)^n h_i s_i + sum_("pair"(i,j)) J_(i,j) s_i s_j \
  ==> quad E(hat(bold(s))) &= -k_b T lr([sum_(i=1)^n h_i s_i + sum_("pair"(i,j)) J_(i,j) s_i s_j]) \
  &= -sum_(i=1)^n (T dot h_i) s_i - sum_("pair"(i,j)) (T dot J_(i,j) ) s_i s_j
  $


  The partition function $Z$ is the normalisation factor for the Boltzmann distribution:

  $
  &sum_(hat(bold(s))) 1/Z exp(-beta E(hat(bold(s)))) = 1 \
  ==> quad &Z = sum_(hat(bold(s))) exp(-beta E(hat(bold(s))))
  $

  For this energy function $Z$ has the form:

  $
  Z = sum_(hat(bold(s))) exp(sum_(i=1)^n h_i s_i + sum_("pair"(i,j)) J_(i,j) s_i s_j)
  $

  #highlight[I started trying to derive a closed-form for this, but was taking a while so will come back to it.]

  ]

+ Take a spin $s_i$: if $h_i$ is positive, which direction will $s_i$ tend to turn to, i.e., which direction of $s_i$ will minimize the associated energy $-h_i s_i$? Take a pair of spins $s_i$ and $s_j$: if $J_(i j)$ is positive, which configurations of $(s_i, s_j)$ minimize the coupling energy $-J_(i j) s_i s_j$?

  Assume that we have inferred the best parameters $h_i$ and $J_(i j)$ for the US supreme court dataset discussed in section 2. How would you interpret the sign of the inferred parameters $h_i$ and $J_(i j)$ in this context?
  #solution[
  Consider a spin $s_i$. If $h_i$ is positive, then the component of the spin's energy attributable to $h_i$ is minimised when $s_i > 0$, such that $-h_i s_i < 0$. Similarly, take two spins $s_i, s_j$ with $i != j$, then if $J_(i j) > 0$ the energy attributable to the spins' interaction is minimised when $"sign"(s_i) = "sign"(s_j)$, such that $s_i s_j > 0$ and $-J_(i j) s_i s_j < 0$. 

  If $h_i$ or $J_(i j)$ were negativethen the opposite result holds ($s_i < 0$, or $"sign"(s_i) != "sign"(s_j)$ respectively). If $h_i = 0$ then $s_i$ has no preferred direction, i.e., the energy due to the field is minimised for any $s_i$. Likewise, if $J_(i j) = 0$ then any configuration of $s_i$ and $s_j$ minimises their interaction energy.

  Suppose now that we have inferred the optimal parameters $h_i$ and $J_(i j)$ for the US supreme court dataset. We can interpret the sign of each parameter by considering its effect _in absence of the other's effect_. For instance, $"sign"(h_i)$ signifies the *political leaning* of the $i$'th judge's votes, in the absence of interactions with other judges. Analogously, $"sign"(J_(i j))$ signifies the *tendency for $i$ and $j$ to vote identically*, in the absence of their individual political leanings. 

  Since we take $+1$ to represent a conservative vote and $-1$ a liberal one, $h_i < 0$ indicates that $i$ has a tendency to vote liberally, and vice versa. Judges $i$ and $j$ tend to vote similarly if $J_(i j) > 0$, and differently if $J_(i j) < 0$. 

  Finally, $h_i = 0$ indicates no particular tendency to vote conservative or liberal, and $J_(i j) = 0$ implies no correlation between $i$ and $j$.
]

== Observables

+ Given a stationary probability distribution of the state $p_g (bold(s))$, what are the definitions of $angle.l s_i angle.r$ and of $angle.l s_i s_j angle.r$?
  #solution[
  For clarity define $s_i$ as a function which extracts the $i$'th element of a vector $s$, that is $s_i: hat(bold(s)) mapsto (hat(bold(s)))_i$. Then $angle.l s_i angle.r$ is the average local magnetisation of the $i$'th spin:

  $
  angle.l s_i angle.r = ensembleavg(s_i)
  $ <eq:avg-local-magnetisation>

  And $angle.l s_i s_j angle.r$ is the average local correlation between the $i$'th and $j$'th spins:

  $
  angle.l s_i s_j angle.r = sum_microstate(s) s_i (microstate(s)) s_j (microstate(s)) dot P(microstate(s))
  $ <eq:avg-local-correlation>
  ]

+ Consider a dataset $hat(bold(s))$ composed of $N$ independent observations of the spins: $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$. Let us denote by $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ the empirical averages of $s_i$ and of $s_i s_j$ respectively (i.e., their average values in the dataset). How would you compute $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ from the data?
  #solution[
    Let $D$ denote the 'dataset', i.e., the set of observation ${microstate(s)^((1)), ..., microstate(s)^((N))}$. We formalise the notion of the data distribution $P_D$ as the proportion of $D$ comprising observations of a particular microstate:
    
    $
    P_D (microstate(s)) = (\# {microstate(s)^((k)) in D | microstate(s)^((k)) = microstate(s)})/N
    $ <eq:data-probability>

    We may define an estimate for $angle.l s_i angle.r$ using $P_D$, such that it can be computed:

    $
    angle.l s_i angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) dot P_D (microstate(s)) \
    &= sum_(microstate(s) in.not D) s_i (microstate(s)) dot 0 + sum_(microstate(s) in D) s_i (microstate(s)) dot P_D (microstate(s)) \
    &= 1/N sum_(microstate(s) in D) s_i (microstate(s)) dot \# {microstate(s)^((k)) in D | microstate(s)^((k)) = microstate(s)} \
    &= 1/N sum_(k=1)^N s_i (microstate(s)^((k)))
    $ <eq:local-magnetisation-using-pdata>


    i.e., the empirical average local magnetisation of spin $i$ can be calculated as the average value of $s_i$ as observed in the data. Note that the sum index, $microstate(s) in D$ on line 3 is taken to mean "$microstate(s)$ occurs in the dataset", rather than being an enumeration of the rows. 

    Likewise for $angle.l s_i s_j angle.r_D$:

    $
    angle.l s_i s_j angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) s_j (microstate(s)) dot P_D (microstate(s)) \
    &= 1/N sum_(microstate(s) in D) s_i (microstate(s)) s_j (microstate(s)) dot \# {microstate(s)^((k)) in D | microstate(s)^((k)) = microstate(s)} \
    &= 1/N sum_(k=1)^N s_i (microstate(s)^((k))) s_j (microstate(s)^((k)))
    $
  ]

+ Assume that the data is stationary and that eah datapoint has been randomly sampled from $p(bold(s))$. Can you show that the empirical averages, $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$, converge to the model averages, respectively $angle.l s_i angle.r$ and $angle.l s_i s_j angle.r$, as the number $N$ of datapoints goes to infinity? (very large dataset)
  #solution[
  In the previous question we defined $P_D (microstate(s))$, for some microstate $microstate(s)$, as the proportion of observations in the dataset $D$ which were equal to $microstate(s)$. We can equivalently define $P_D (microstate(s))$ as an average over the function $bb(1)[microstate(s)^((k)) = microstate(s)]$, which is 1 if this condition holds, and 0 otherwise:
  $
  P_D (microstate(s)) = 1/N sum_(k=1)^N bb(1)[microstate(s)^((k)) = microstate(s)]
  $
  By the Law of Large Numbers, as $n -> infinity$, this average converges to its expected value:
  $
  lim_(N -> infinity) P_D (microstate(s)) = P(microstate(s))
  $
  The main result follows naturally by considering $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ as $N -> infinity$:
  $
  lim_(N -> infinity) angle.l s_i angle.r_D &= lim_(N -> infinity) lr([sum_(microstate(s) in D) s_i (microstate(s)) dot P_D (microstate(s))]) 
  &= sum_(microstate(s)) s_i (microstate(s)) dot P(microstate(s)) 
  &= angle.l s_i angle.r
    $
  And for the local correlation:
  $
  lim_(N -> infinity) angle.l s_i s_j angle.r_D &= lim_(N -> infinity) lr([sum_(microstate(s) in D) s_i (microstate(s)) s_j (microstate(s)) dot P_D (microstate(s))])
  &= sum_(microstate(s)) s_i (microstate(s)) s_j (microstate(s)) dot P(microstate(s))
  &= angle.l s_i s_j angle.r
  $
  ]

== Maximum Entropy models

+ Consider a spin system with stationary probability distribution $p(bold(s))$. Can you recall the definition of the Shannon entropy $S[p(bold(s))]$? As mentioned above for the Boltzmann distribution, we will take $k_b = 1$.
  
  #solution[
  With $k_b = 1$, the Shannon entropy is:
  $
  S[p(bold(s))] = -sum_(microstate(s)) p(microstate(s)) log p(microstate(s))
  $
  ]

  The Ising model in Eq. (1) can be seen as a _Maximum Entropy Model_, constrained to reproduce the data local magnetisation and local correlation, i.e., constrained to reproduce all the data averages $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ (for all spins $s_i$ and $s_j$). We also want $p(bold(s))$ to be normalised, which introduces the additional constraint $sum_bold(s) p(bold(s)) = 1$. To summarise, we are looking for the set of $2^n$ probabilities $p(bold(s))$ such that $S[p(bold(s))]$ is maximal, and such that 
  $
  sum_bold(s) p(bold(s)) = 1 quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) = angle.l s_i angle.r_D quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) = angle.l s_i s_j angle.r_D
  $
  where $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ are constants that are computed from the data for all distinct $s_i$ and $s_j$. Note that to be more precise, we wrote $s_i (bold(s))$ (instead of just $s_i$) to specify that this is the value of $s_i$ in the state $bold(s)$ (this will help with the next questions).

+ How many constraints are there in total?
  #solution[
  The constraints are as follows:
  - *Maximisation:* 1 constraint,
  - *Normalisation:* 1 constraint,
  - *Average local magnetisation:* One per spin ($n$ total)
  - *Average local correlation:* One for each pair of distinct spins ($(n(n-1))/2$ total)

  Thus the total number of constraints is
  $
  1 + 1 + n + (n(n-1))/2 = (n(n + 1) + 1)/2
  $
  ]

  To find the shape of the distributions $p(bold(s))$ that maximises the entropy while satisfying these constraints, we introduce an auxiliary function:
  $
  U[p(bold(s))] = S[p(bold(s))] + lambda_0 (sum_s p(bold(s)) - 1) + sum_(i=1)^n alpha_i (sum_bold(s) p(bold(s)) s_i (bold(s)) - angle.l s_i angle.r_D) \
  + sum_("pair"(i,j))^n eta_(i j) (sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) - angle.l s_i s_j angle.r_D)
  $
  where we have introduced a parameter in front of each constraint we want to impose. These parameters ($lambda_0$, $alpha_i$, and $eta_(i j)$) are called Lagrange multipliers. To find $p(bold(s))$ one must maximise this auxiliary function with respect to the $2^n$ probabilities $p(bold(s))$.

+ Let us fix a choice of a state $bold(s)$. The probability $p_bold(s) = p(bold(s))$ is a parameter of $U[bold(p)]$ where $bold(p)$ is the vector of the $2^n$ probabilities. Can you show that:
  $
  (diff U[bold(p)])/(diff p_s) = -log(p_s) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)) 
  $
  #solution[
  For clarity, we treat the terms in the derivative one at a time. Observe that since we are taking the partial derivative of $U$ with respect to a single element in the vector $p(bold(s))$, in each of the sums over $microstate(s)$ in $U$, all terms will be annhialated by the derivative, with the exception of the particular $bold(s)$ which we have fixed.  

  First, examining the Shannon entropy:
  $
  diff/(diff p_s) S[p(s)] &= diff/(diff p_s) lr([-sum_microstate(s) p(microstate(s)) log p(microstate(s))]) \
  &= -log(p_s) - p_s dot 1/p_s \
  &= -log(p_s) - 1
  $

  Next, the normalisation constraint:
  $
  diff/(diff p_s) lambda_0 lr((sum_microstate(s) p(microstate(s)) - 1)) &= diff/(diff p_s) lambda_0 p_s &= lambda_0
  $

  The local average magnetisation constraint:
  $
  diff/(diff p_s) sum_(i=1)^n alpha_i lr((sum_microstate(s) p(microstate(s)) s_i (microstate(s)) - angle.l s_i angle.r_D)) &= sum_(i=1)^n alpha_i lr((diff/(diff p_s) p_s s_i (bold(s)))) \
  &= sum_(i=1)^n alpha_i s_i (bold(s))
  $

  And lastly the average local correlation constraint:
  $
  diff/(diff p_s) sum_("pair"(i,j)) eta_(i j) lr((sum_microstate(s) p(microstate(s)) s_i (microstate(s)) s_j microstate(s) - angle.l s_i s_j angle.r_D)) &= sum_("pair"(i,j)) eta_(i j) lr((diff/(diff p_s) p_s s_i (bold(s)) s_j (bold(s)))) \
  &= sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))
  $

  Taking these results together, we arrive at the desired result:
  $
  (diff U[bold(p)])/(diff p_s) = -log(p_s) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))
  $
  ]

+ Can you show that the most general expression of $p_bold(s)$ with maximal entropy that satisfying the constraints in Eq. (2) is Eq. (1)? Give the relation between $lambda_0$ and the partition function $Z$. How are the parameters $alpha_i$ or $eta_(i j)$ related to the parameters $h_i$ and $J_(i j)$?
  #solution[
  The constrained optimisation problem described by $U$ is optimised iff, for all microstates $bold(s)$ and spins $i, j$, $ (diff U)/(diff p_s) = (diff U)/(diff lambda_0) = (diff U)/(diff alpha_i) = (diff U)/(diff eta_(i j)) = 0 $
  Note that the partial derivatives with respect to the constraints are zero exactly when those constraints are satisfied. The partial derivative with respect to $p_s$ is zero for a particular $bold(s)$ when
  $
  0 &= -log(p_s) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)) \
  p_s &= exp(-(1 - lambda_0) + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s))) \
  &= 1/exp(1 - lambda_0) exp(sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)))
  $

  Comparing this expression to Eq. (1), we obtain the following equalities:
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

+ Can you show that in that case the model is assuming the variables are independent from each other, i.e., that we can write the joint probability distribution as a product of a probability distribution over each variable $p_bold(g) (bold(s)) = product_(i=1)^n p_(bold(h_i)) (s_i)$? What is the probability distribution $p_(bold(h_i)) (s_i)$ for the spin variable $s_i$?
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

+ Take one of the spin variables $s_i$. We recall that $angle.l s_i angle.r_D$ is the average value of $s_i$ in the data (given a dataset, this quantity is a constant), and that $angle.l s_i angle.r = sum_bold(s) p(bold(s)) s_i$ is the model average of $s_i$. Can you show that the value of the parameter $h_i$ that satisfies the constraint $angle.l s_i angle.r = angle.l s_i angle.r_D$ is:
  $
  h_i = tanh^(-1) (angle.l s_i angle.r_D),
  $
  where $tanh^(-1) (x)$ denotes the inverse of the hyperbolic tangent? In particular, in that case the probability distribution over $s_i$ in the model is exactly equal to the empirical distribution of $s_i$.
  #solution[
  Let $i$ be a spin in a model with $n$ spins, and $angle.l s_i angle.r_D$ the average value of $s_i$ in a dataset $D$. Suppose that $angle.l s_i angle.r = angle.l s_i angle.r_D$ in the model, then from the definition of $angle.l s_i angle.r$
  $
  angle.l s_i angle.r = sum_(microstate(s)) s_i (microstate(s)) dot p_bold(g) (microstate(s)) = angle.l s_i angle.r_D
  $ <eq:2-4-2-definition>
  From the previous exercise, we have that $p_bold(g) (microstate(s)) = inline(product_(j=1)^n) p_bold(h_j) (s_j)$. Taking $hat(bold(s))_i = (s_1, ..., s_(i-1), s_(i+1), ..., s_n)$ to denote the vector of spin variables excluding $s_i$, we use this result to rewrite @eq:2-4-2-definition in terms of $p_bold(h_i)$:
  $
  angle.l s_i angle.r_D &= sum_(microstate(s)) s_i (microstate(s)) dot p_bold(h_i) (s_i (microstate(s))) dot product_(j != i) p_bold(h_j) (s_j (microstate(s))) \
  &= sum_(s in plus.minus 1) s p_bold(h_i) (s) dot overbrace(sum_(hat(bold(s))_(-i)) product_(j != i) p_bold(h_j (s_j (hat(bold(s))_(-i)))), = 1) \
  &= p_bold(h_i) (1) - p_bold(h_i) (-1) \
  &= (exp(h_i) - exp(-h_i))/(2cosh(h_i)) \
  &= sinh(h_i)/cosh(h_i) \
  &= tanh(h_i)
  $

  In the second line of reasoning we group the original summation terms by their $s_i$ component. The inner summation is then the sum over the probabilities of observing any other microstate for the remaining $n-1$ variables. The normalisation condition of probability distributions implies that this is equal to 1.

  Finally, we obtain the desired result, $h_i = tanh^(-1)(angle.l s_i angle.r_D)$. 
  ]

+ In Eq. (6), we observe that:
  - If $angle.l s_i angle.r_D > 0$, then the inferred $h_i$ is also positive;
  - Reciprocally, $angle.l s_i angle.r_D < 0$, then the inferred $h_i$ is also negative.
  How does this connect with the tendency of the $i$'th judge to vote on average more liberal or more conservative? Is this result coherent with the general comments that we did in Question Q1.4.?
  #solution[
  This result is coherent with our comments in Question Q1.4, in which we interpreted the sign of $h_i$ as reflecting the sign of $i$'s political leaning ($-1$ for liberal, $+1$ for conservative). As $angle.l s_i angle.r_D$ is the average value of $s_i$ in a dataset $D$, a positive value occurs when $i$ votes conservative in more than 50% of cases. Likewise, a negative value occurs when $i$ votes liberal in more than 50% of cases. 

  As in Q1.4, if the $i$'th judge has an equal number of liberal and conservative votes in $D$, then $angle.l s_i angle.r_D = h_i = 0$. 

  While in Q1.4 our discussion was concerned with a model which included interactions, the comments are still relevant here, as we interpreted the sign of $h_i$ to reflect political leaning in absence of interactions with other judges. In this model we simply make this assumption explicit by taking $J_(i j) = 0$.
  ]

== Statistical inference: maximising the log-likelihood function

*Introducing the likelihood function.* Looking more closely at Eq. (1), one can see that it does not just define a single probability distribution, but many of them: there is one probability distribution for each value of the set of parameters $bold(g)$. More precisely, the distribution in Eq. (1) changes continuously as one continuously varies the parameters in $bold(g)$. We say that Eq. (1) defines a _parametric family of probability distributions_. The inference procedure consists in finding the value of the parameters $bold(g)$ that maximises the probability that the model $p_bold(g) (bold(s))$ produces the data.

To do so, we introduce the _log-likelihood function_:
$
cal(L)(bold(g)) = log P_bold(g) (microstate(s))
$ <eq:2-5-log-likelihood>
where $P_bold(g) (microstate(s))$ is the probability that the model $p_bold(g) (bold(s))$ produces the dataset $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$. Note that $cal(L)(bold(g))$ is a function of the parameters $bold(g)$. The inference procedure therefore consists in finding the value $bold(g)^star$ of the parameters that maximises $cal(L)(bold(g))$. For the moment we will assume that there exists only a unique such value of $bold(g)$.

+ We assume that, in the dataset $hat(bold(s))$, all datapoints are independently samples from an underlying distribution $p_bold(g) (bold(s))$. Can you show that the log-likelihood function can be re-written as:
  $ cal(L) = N sum_(bold(s)) p_D (bold(s)) log p_bold(g) (bold(s))$
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
  ]

  *Ising model.* We now take the model distribution $p_bold(g) (bold(s))$ to be given by the Ising model in Eq. (1).
+ Taking the first derivative of $cal(L)(bold(g))$ with respect to a parameter $h_i$, can you show that at the maximum of $cal(L)(bold(g))$ we have that $angle.l s_i angle.r = angle.l s_i angle.r_D$? Similarly, taking the first derivative of $cal(L)(bold(g))$ with respect to a parameter $J_(i j)$, can you show that at the maximum of $cal(L)(bold(g))$ we have that $angle.l s_i s_j angle.r = angle.l s_i s_j angle.r_D$?
  #solution[
  Let $p_bold(g) (bold(s))$ be as defined in Eq. (1), and let $i$ correspond to the spin $s_i$. Then $cal(L)(bold(g))$ is:
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

+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]
+
  #solution[]

