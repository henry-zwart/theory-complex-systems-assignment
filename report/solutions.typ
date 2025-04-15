#import "@preview/cetz:0.3.4"
#import "@preview/fletcher:0.5.7" as fletcher: diagram, node, edge, shapes
#import "@local/uva-game-theory:0.1.0": *


= Modelling the activity of a single neuron

#let neuron_results = json("../results/neuron_modelling.json")

+ Can you plot the distribution $P(tau)$ of the time intervals $tau$ between successive spikes? Check that there is indeed a rafactory period, i.e., a time interval $tau_0$ after each spike, during which the neuron doesn't spike again. What is the duration $tau_0$ for this time interval?
  #solution[
  #figure(image("../results/figures/neuron_interspiking_distribution.svg"))

  Refractory period duration $tau_0 = #neuron_results.tau0 m s$ calculated as the minimum observed duration between neuron spikes.
  ]

+ Can you check that the decay of the distribution $P(tau)$ of inter-spike intervals is indeed exponential? Measure the corresponding decay rate $lambda$
  #solution[]

+ Can you deduce an analytical expression for the distribution of inter-spike time interval $P(tau)$ of the delayed Poisson process as a function of $lambda$ and $tau_0$? Compare your model distribution to the one obtained from the data.
  #solution[]

+ Using your model, can you generate another 1000 (spike times) datapoints?
  #solution[]

+ What is the average spiking rate $f$ of the neuron in the data? How is $f$ analytically related to $tau_0$ and $lambda$ that you have previously measured?
  #solution[]

= Modelling binary data with the Ising model

== Pairwise spin model

+ How many terms are in the sum over the $"pair"(i, j)$? Can you deduce what is the number of parameters in the vector $g = (h_1,..., h_n, J_(1,2),..., J_n$? Can you re-write the sum over the $"pair"(i, j)$ as a double sum over $i$ and $j$ (without counting twice each pair)?
  #solution[]

+ Can you write down explicitly the terms in the exponential of Eq. (1) for a system with $n = 3$ spins?
  #solution[]

+ In Eq. (1), we can recognize the Boltzmann distribution, in which the parameter $beta = 1\/(k_b T)$ was taken equal to 1 (more precisely, the constant $k_B$ was taken equal to 1, and the temperature parameter $T$ was absorbed in the parameters $h_i$ and $J_(i j)$). What is the energy function associated with the Boltzmann distribution in that case? What is the partition function and what is its general expression?
  #solution[]

+ Take a spin $s_i$: if $h_i$ is positive, which direction will $s_i$ tend to turn to, i.e., which direction of $s_i$ will minimize the associated energy $-h_i s_i$? Take a pair of spins $s_i$ and $s_j$: if $J_(i j)$ is positive, which configurations of $(s_i, s_j)$ minimize the coupling energy $-J_(i j) s_i s_j$?

  Assume that we have inferred the best parameters $h_i$ and $J_(i j)$ for the US supreme court dataset discussed in section 2. How would you interpret the sign of the inferred parameters $h_i$ and $J_(i j)$ in this context?
  #solution[]

== Observables

+ Given a stationary probability distribution of the state $p_g(bold(s))$, what are the definitions of $angle.l s_i angle.r$ and of $angle.l s_i s_j angle.r$?
  #solution[]

+ Consider a dataset $hat(bold(s))$ composed of $N$ independent observations of the spins: $hat(bold(s)) = (bold(s)^((1)), ..., bold(s)^((N)))$. Let us denote by $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ the empirical averages of $s_i$ and of $s_i s_j$ respectively (i.e., their average values in the dataset). How would you compute $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ from the data?
  #solution[]

+ Assume that the data is stationary and that eah datapoint has been randomly sampled from $p(bold(s))$. Can you show that the empirical averages, $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$, converge to the model averages, respectively $angle.l s_i angle.r$ and $angle.l s_i s_j angle.r$, as the number $N$ of datapoints goes to infinity? (very large dataset)
  #solution[]

== Maximum Entropy models

+ Consider a spin system with stationary probability distribution $p(bold(s))$. Can you recall the definition of the Shannon entropy $S[p(bold(s))]$? As mentioned above for the Boltzmann distribution, we will take $k_b = 1$.
  
  The Ising model in Eq. (1) can be seen as a _Maximum Entropy Model_, constrained to reproduce the data local magnetisation and local correlation, i.e., constrained to reproduce all the data averages $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ (for all spins $s_i$ and $s_j$). We also want $p(bold(s))$ to be normalised, which introduces the additional constraint $sum_bold(s) p(bold(s)) = 1$. To summarise, we are looking for the set of $2^n$ probabilities $p(bold(s))$ such that $S[p(bold(s))]$ is maximal, and such that 
  $
  sum_bold(s) p(bold(s)) = 1 quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) = angle.l s_i angle.r_D quad "and" quad sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) = angle.l s_i s_j angle.r_D
  $
  where $angle.l s_i angle.r_D$ and $angle.l s_i s_j angle.r_D$ are constants that are computed from the data for all distinct $s_i$ and $s_j$. Note that to be more precise, we wrote $s_i (bold(s))$ (instead of just $s_i$) to specify that this is the value of $s_i$ in the state $bold(s)$ (this will help with the next questions).
  #solution[]

+ How many constraints are there in total?

  To find the shape of the distributions $p(bold(s))$ that maximises the entropy while satisfying these constraints, we introduce an auxiliary function:
  $
  U[p(bold(s))] = S[p(bold(s))] + lambda_0 (sum_s p(bold(s)) - 1) + sum_(i=1)^n alpha_i (sum_bold(s) p(bold(s)) s_i (bold(s)) - angle.l s_i angle.r_D) \
  + sum_("pair"(i,j))^n eta_(i j) (sum_bold(s) p(bold(s)) s_i (bold(s)) s_j (bold(s)) - angle.l s_i s_j angle.r_D)
  $
  where we have introduced a parameter in front of each constraint we want to impose. These parameters ($lambda_0$, $alpha_i$, and $eta_(i j)$) are called Lagrange multipliers. To find $p(bold(s))$ one must maximise this auxiliary function with respect to the $2^n$ probabilities $p(bold(s))$.
  #solution[]

+ Let us fix a choice of a state $bold(s)$. The probability $p_bold(s) = p(bold(s))$ is a parameter of $U[bold(p)]$ where $bold(p)$ is the vector of the $2^n$ probabilities. Can you show that:
  $
  (diff U[bold(p)])/(diff p_s) = -log(p_s) - 1 + lambda_0 + sum_(i=1)^n alpha_i s_i (bold(s)) + sum_("pair"(i,j)) eta_(i j) s_i (bold(s)) s_j (bold(s)) 
  $
  #solution[]

+ Can you show that the most general expression of $p_bold(s)$ with maximal entropy that satisfying the constraints in Eq. (2) is Eq. (1)? Give the relation between $lambda_0$ and the partition function $Z$. How are the parameters $alpha_i$ or $eta_(i j)$ related to the parameters $h_i$ and $J_(i j)$?
  #solution[]

== Statistical inference: model with no couplings

+
  #solution[]
+
  #solution[]
+
  #solution[]

== Statistical inference: maximising the log-likelihood function

+
  #solution[]
+
  #solution[]

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

