import pyprob
from pyprob import Model
from pyprob.distributions import Normal

import torch
import numpy as np
import math
import matplotlib.pyplot as plt

class GaussianUnknownMean(Model):
    def __init__(self):
        super().__init__(name='Gaussian with unknown mean') # give the model a name
        self.prior_mean = 1
        self.prior_std = math.sqrt(5)
        self.likelihood_std = math.sqrt(2)

    def forward(self): # Needed to specifcy how the generative model is run forward
        # sample the (latent) mean variable to be inferred:
        mu = pyprob.sample(Normal(self.prior_mean, self.prior_std)) # NOTE: sample -> denotes latent variables

        # define the likelihood
        likelihood = Normal(mu, self.likelihood_std)

        # Lets add two observed variables
        # -> the 'name' argument is used later to assignment values:
        pyprob.observe(likelihood, name='obs0') # NOTE: observe -> denotes observable variables
        pyprob.observe(likelihood, name='obs1')

        # return the latent quantity of interest
        return mu

model = GaussianUnknownMean()
model.learn_inference_network(num_traces=20000,
        observe_embeddings={'obs0' : {'dim' : 32},
            'obs1': {'dim' : 32}},
        inference_network=pyprob.InferenceNetwork.LSTM)

observed_list = [12, 10] # Observations
# sample from posterior (200 samples)
posterior = model.posterior_results(
        num_traces=200, # the number of samples estimating the posterior
        inference_engine=pyprob.InferenceEngine.IMPORTANCE_SAMPLING_WITH_INFERENCE_NETWORK, # specify which inference engine to use
        observe={
            'obs0': observed_list[0],
            'obs1': observed_list[1]} # assign values to the observed values
        )
print(posterior.expectation(lambda x: x))
