from abc import abstractmethod
from typing import Union, List

import torch
import torch.nn as nn
import torch.nn.functional as F


class TinyVAE(nn.Module):
    def __init__(self, encoder: nn.Module = None, decoder: nn.Module = None):
        super().__init__()

        self.encoder = encoder
        self.decoder = decoder

    def encode(self, x):
        return self.encoder.forward(x)

    def decode(self, latent_x):
        return self.decoder.forward(latent_x)

    def forward(self, x):
        latent_x, mu, logvar = self.encoder.forward(x)  # encoding
        recon_x = self.decoder.forward(latent_x)  # decoding

        return recon_x, mu, logvar


class LinearModel(nn.Module):
    def __init__(
        self,
        input_size: int = 50,
        output_size: int = 10,
        num_hidden_layers: int = 2,
        hidden_layer_size: int = 100,
        dropout: Union[float, List[float]] = 0.85,
    ):
        """
        Initialize a `num_hidden_layers` + input layer + outputs layer multilayer perceptron for encoder.

        Parameters
        ----------
        input_size: int
            Length of linear numeric input vector (input is embedding for encoders and latent space vector for decoder)

        output_size: int
            Length of linear numeric output vector (output is latent space for encoders and embedding vector for decoder)

        num_hidden_layers: int
            Number of equivalent hidden layers

        hidden_layer_size: int
            Number of neurons in each hidden layer

        droupout: float or list of floats
            Dropout rate (d) of the layers, (no dropout) 0 ≤ d < 1 (all dropout)
        """
        super().__init__()

        self.in_sz = input_size
        self.out_sz = output_size
        self.num_hid = num_hidden_layers
        self.hid_sz = hidden_layer_size
        self.num_layers = self.num_hid + 2
        if type(dropout) == float:
            self.dropout_p = [dropout] * (self.hid_sz + 1)
        elif len(dropout) == self.hid_sz + 1:
            self.dropout_p = dropout
        else:
            raise ValueError(
                f"Input `dropout` is of size {len(dropout)} but has to be {self.hid_sz + 1}"
            )

        self.layers = nn.ModuleList()  # initialize list of layers
        self.dropout = nn.ModuleList()  # initialize list of dropout-able layers

        # stacking neural network layers
        self.layers.append(nn.Linear(self.in_sz, self.hid_sz))
        self.dropout.append(nn.Dropout(p=self.dropout_p[0]))  # dropout for input layer

        for i in range(self.num_hid):
            self.dropout.append(
                nn.Dropout(p=self.dropout_p[i])
            )  # dropout for hidden layers
            self.layers.append(
                nn.Linear(self.hid_sz, self.hid_sz)
            )  # fully-connected hidden layers

        self.layers.append(nn.Linear(self.hid_sz, self.out_sz))  # output layer

    @abstractmethod
    def forward(self, x):
        pass


class LinearEncoder(LinearModel):
    def __init__(
        self,
        embed_size: int = 50,
        latent_size: int = 10,
        num_hidden_layers: int = 2,
        hidden_layer_size: int = 100,
        dropout: Union[float, List[float]] = 0.85,
    ):
        super().__init__(
            input_size=embed_size,
            output_size=latent_size,
            num_hidden_layers=num_hidden_layers,
            hidden_layer_size=hidden_layer_size,
            dropout=dropout,
        )

    def reparam(self, mu, log_var):
        """
        Parameters
        ----------
        Use reparameterization trick to allow gradient flow

        mu: int
            Mean in encoder's latent space

        log_var:
            Log of variance in enconder's latent space
        """
        std = torch.exp(0.5 * log_var)  # standard deviation
        eps = torch.randn_like(std)  # `randn_like` as we need the same size
        sample = mu + (eps * std)  # sampling as if coming from the input space

        return sample

    def forward(self, x):
        for i in range(self.num_hid + 1):  # excluding output layer
            x = F.relu(self.dropout[i](self.layers[i](x)))
        x = self.layers[-1](x).view(-1, 2, self.out_sz // 2)

        mu = x[:, 0, :]  # the first feature values as mean
        log_var = x[:, 1, :]  # the other feature values as variance
        latent_x = self.reparam(mu, log_var)

        return latent_x, mu, log_var


class LinearDecoder(LinearModel):
    def __init__(
        self,
        embed_size: int = 50,
        latent_size: int = 10,
        num_hidden_layers: int = 2,
        hidden_layer_size: int = 100,
        dropout: Union[float, List[float]] = 0.85,
    ):
        """
        Initialize a `num_hidden_layers` + input layer + outputs layer multilayer perceptron for decoder.

        Parameters
        ----------
        input_size: int
            Length of input vector, should be a linear numeric vector

        num_hidden_layers: int
            Number of equivalent hidden layers

        hidden_layer_size: int
            Number of neurons in each hidden layer

        droupout: float or list of floats
            Dropout rate (d) of the layers, (no dropout) 0 ≤ d < 1 (all dropout)
        """
        super().__init__(
            input_size=latent_size,
            output_size=embed_size,
            num_hidden_layers=num_hidden_layers,
            hidden_layer_size=hidden_layer_size,
            dropout=dropout,
        )

    def forward(self, latent_x):
        x = latent_x
        for i in range(self.num_hid + 1):  # excluding output layer
            x = F.relu(self.dropout[i](self.layers[i](x)))
        recon_x = torch.sigmoid(self.layers[-1](x))

        return recon_x


class ConvolutionalVAE(nn.Module):
    def __init__(
        self,
        input_size=308,
        num_hidden_layers=2,
        hidden_layer_size=100,
        dropout_input=0.85,
        dropout_hidden=0.65,
        conv=False,
        num_conv_layers=1,
        conv_filt_sz=5,
        conv_stride=2,
    ):
        pass