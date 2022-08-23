import pyprob
from pyprob import Model
from pyprob.distributions import Normal

import numpy as np
import torch
from pyprob.distributions import Gamma, Mixture

class GMM(Model):
    def __init__(self, K=2):
        super().__init__(name='Gaussian mixture model')
        self.d = 2
        self.K = K

    def forward(self):
        muxs = list()
        muys = list()
        probs = [1./self.K for _ in range(self.K)]
        for _ in range(self.K):
            muxs.append(pyprob.sample(Normal(0, 1)))
            muys.append(pyprob.sample(Normal(0, 1)))

        # define the likelihood factorized over coordinates
        likelihoodx = Mixture(
                [Normal(muxs[i], 0.1) for i in range(self.K)],
                probs=probs)
        likelihoody = Mixture(
                [Normal(muys[i], 0.1) for i in range(self.K)],
                probs=probs)

        pyprob.observe(likelihoodx, name='obs0x')
        pyprob.observe(likelihoody, name='obs0y')
        pyprob.observe(likelihoodx, name='obs1x')
        pyprob.observe(likelihoody, name='obs1y')

        return torch.tensor([muxs, muys]).T


if __name__ == '__main__':
    model = GMM()
    # prior = model.prior_results(num_traces=10)
    # print(prior.get_values())
    model.learn_inference_network(
            num_traces=20000,
            observe_embeddings={
                'obs0x' : {'dim' : 16},
                'obs0y' : {'dim' : 16},
                'obs1x' : {'dim' : 16},
                'obs1y' : {'dim' : 16},
            },
            inference_network=pyprob.InferenceNetwork.LSTM)

    observed_list = [100, 100, -100, -100] # Observations
    # sample from posterior (200 samples)
    posterior = model.posterior_results(
            num_traces=200, # the number of samples estimating the posterior
            inference_engine=pyprob.InferenceEngine.IMPORTANCE_SAMPLING_WITH_INFERENCE_NETWORK, # specify which inference engine to use
            observe={
                'obs0x': observed_list[0],
                'obs0y': observed_list[1],
                'obs1x': observed_list[2],
                'obs1y': observed_list[3]} # assign values to the observed values
            )
    print(posterior.expectation(lambda x: x))
