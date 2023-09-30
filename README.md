# Estimating Dynamic Treatment Regimes through Deep Reinforcement Learning in simulated cancer trial scenarios

## **Description**
This repo contains an implementation of Deep Q-Network and Deep Recurrent Q-Network (DQN) considering different models (RNN, LSTM and DNN) in a Dynamic Treatment Regime scenario (DTR). Physicians must be able to make numerous healthcare decisions during the course of a patient’s illness by monitoring both the individual’s medical history and the evolution of the disease together with its possible consequences. DTR is a treatment design where a set of decision rules aim to find the optimal regime that, if followed by the individual, would yield the best results on average. DQN can be a viable data driven solution to provide support in this scheme.

## **Approach**
The approach chosen to evaluate DRL techniques in DTR contexts is a Deep Q-Network using three type of Neural Network: RNN, LSTM and DNN. The RL technique used as a benchmark is instead QL via OLS. Several strategies have been used for DRL methods, which are discussed in more detail in the papers _Exploration Policies_ and _Target Network Updates_. Every approach can be enhanced with several exploration strategies, like deterministic epsilon-greedy and softmax.

## **Simulative Models**
Given the limited availability of longitudinal observational data describing clinical trials due to privacy constraints and contractual regulations, a common way to obtain this data is to resort to simulation. In particular, this work recall Kyle Humphrey’s simulation study which in turn was largely inspired by Yufan Zhao, Kosorok, and Zeng et al. and Yufan Zhao, Zeng et al. They describe a generic cancer clinical trial for a predefined number of stages with the aim of finding the optimal treatment regime (dose of a drug) for each patient maximizing survival time.

## **Papers**
In the _Papers_ folder you'll find a collection of detailed PDFs that delve into the intricacies of the DRL approaches I've employed, complete with thorough references to the bibliography.
