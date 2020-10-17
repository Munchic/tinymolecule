from typing import Union, List

import torch
import torch.nn as nn
import torch.nn.functional as F

class LinearVAE(nn.Module):
    def __init__(
        self,
        encoder: nn.Module = LinearEncoder
        decoder: nn.Module = LinearDecoder
    ):
        super(LinearVAE, self).__init__()

        self.inp_sz = input_size
        self.num_hid = num_hidden_layers
        self.hid_sz = hidden_layer_size
        self.num_conv_layers = num_conv_layers
        self.conv_filt_sz = conv_filt_sz
        self.conv_stride = conv_stride
        self.hidden = nn.ModuleList()  # initialize list of layers
        self.conv = nn.ModuleList()
        self.dropout = nn.ModuleList()  # initialize list of dropout-able layers
        self.conv_flag = conv

        # encoder
        self.enc1 = nn.Linear(in_features=784, out_features=512)
        self.enc2 = nn.Linear(in_features=512, out_features=features * 2)

        # decoder
        self.dec1 = nn.Linear(in_features=features, out_features=512)
        self.dec2 = nn.Linear(in_features=512, out_features=784)

    def reparameterize(self, mu, log_var):
        """
        :param mu: mean from the encoder's latent space
        :param log_var: log variance from the encoder's latent space
        """
        std = torch.exp(0.5 * log_var)  # standard deviation
        eps = torch.randn_like(std)  # `randn_like` as we need the same size
        sample = mu + (eps * std)  # sampling as if coming from the input space

        return sample

    def forward(self, x):
        # encoding
        x = F.relu(self.enc1(x))
        x = self.enc2(x).view(-1, 2, features)
        # get `mu` and `log_var`
        mu = x[:, 0, :]  # the first feature values as mean
        log_var = x[:, 1, :]  # the other feature values as variance
        # get the latent vector through reparameterization
        z = self.reparameterize(mu, log_var)

        # decoding


        return reconstruction, mu, log_var


class LinearModel(nn.Module):
    def __init__(
        self,
        input_size: int = 50,
        output_size: int = 10
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
        super(LinearModel, self).__init__()

        self.in_sz = input_size
        self.out_sz = input_size
        self.num_hid = num_hidden_layers
        self.hid_sz = hidden_layer_size
        self.num_layers = self.num_hid + 2
        if type(dropout) == float:
            self.dropout_p = [dropout] * (self.hid_sz + 1)
        elif len(dropout) == self.hid_sz + 1:
            self.dropout_p = dropout
        else:
            raise ValueError(f"Input `dropout` is of size {len(dropout)} but has to be {self.hid_sz + 1}")
        
        self.layers = nn.ModuleList()  # initialize list of layers
        self.dropout = nn.ModuleList()  # initialize list of dropout-able layers

        # stacking neural network layers
        self.dropout.append(nn.Dropout(p=dropout_p[0]))  # dropout for input layer
        self.layers.append(
            nn.Linear(self.in_sz, self.hid_sz)
        )

        for i in range(1, self.num_hid):
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
        latent_size: int = 10
        num_hidden_layers: int = 2,
        hidden_layer_size: int = 100,
        dropout: Union[float, List[float]] = 0.85,
    ):
        super(LinearEncoder, self).__init__(
            input_size=embed_size,
            output_size=latent_size,
            num_hidden_layers=num_hidden_layers,
            hidden_layer_size=hidden_layer_size,
            dropout=dropout
        )

    def reparam(self, mu, log_var)
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
        for i in range(self.num_layers - 1):  # excluding output layer
            x = F.relu(self.dropout[i](self.layers[i](x)))
        x = self.layers[-1](x).view(-1, 2, self.output_size)
        
        mu = x[:, 0, :]  # the first feature values as mean
        log_var = x[:, 1, :]  # the other feature values as variance
        latent_x = self.reparam(mu, log_var)

        return latent_x

class LinearDecoder(nn.Module):
    def __init__(
        self,
        input_size: int = 50,
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

    def forward(self, latent_x):
        x = latent_x
        for i in range(self.num_layers - 1):  # excluding output layer
            x = F.relu(self.dropout[i](self.layers[i](x)))
        x = F.relu(self.dec1(latent_x))
        reconstruction = torch.sigmoid(self.dec2(x))


class ConvolutionalVAE(LinearVAE):
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