import os
import pickle
from itertools import permutations, product
import copy
import time

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.style.use("ggplot")


def simple_train(
    model, molecloader, criterion, optimizer, num_epochs=5, show_output=True
):

    since = time.time()

    for e in range(num_epochs):
        if show_output:
            print(
                "epoch {}/{} started at {:.4f} s".format(
                    e + 1, num_epochs, time.time() - since
                )
            )

        model.train()
        running_loss = 0.0

        for _molec in molecloader:
            molec = _molec["tokenized"]
            optimizer.zero_grad()

            with torch.set_grad_enabled(True):
                molec_recon, mu, logvar = model(molec)

                ce_loss = criterion(F.normalize(molec_recon), F.normalize(molec))
                loss = vae_loss(ce_loss, mu, logvar)
                running_loss += loss.item()

                loss.backward()
                optimizer.step()

        train_loss = running_loss / len(molecloader)
        print("epoch_loss:", train_loss)


def vae_loss(ce_loss, mu, logvar):
    # full reconstruction loss
    ce = ce_loss
    kd = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return ce + kd


def trainVAE(
    model,
    molecloader,
    criterion,
    optimizer,
    scheduler=None,
    num_epochs=25,
    learning_rate=1e-2,
    show_output=True,
):
    """
    Train neural networks in VAE using gradient descent and backprop

    Parameters
    ----------
    model: nn.Module
        Defined neural network
    molecloader: dict
        Dictionary of torch.utils.data.DataLoader for train, val

    Returns
    ----------
    model: nn.Module
        Trained neural network
    train_history: dict
        Dictionary of defined metrics for train and val sets for each epoch
    state_dicts: list
        List of model state dictionaries at each epoch
    """

    since = time.time()
    train_history = init_train_hist()
    best_model_wts = copy.deepcopy(model.state_dict())
    best_auc = 0.0
    state_dicts = []

    phases = ["train"]

    for e in range(num_epochs):
        if show_output:
            print(
                "epoch {}/{} started at {:.4f} s".format(
                    e + 1, num_epochs, time.time() - since
                )
            )
        for phase in phases:
            if phase == "train":
                model.train()  # set model to training mode
            else:
                model.eval()  # set model to evaluation mode

            count = 0  # just use length of data
            running_loss = 0.0

            # storage for targets and predictions from batches
            y_actual = []
            y_pred = []
            y_proba = []

            for item in molecloader[phase]:
                pep = item["peptide"]
                target = item["target"].view(-1, 1)
                optimizer.zero_grad()

                with torch.set_grad_enabled(phase == "train"):
                    output = model(pep.float())
                    _pos = torch.ones(output.size())
                    _neg = torch.zeros(output.size())
                    preds = torch.where(output > 0.5, _pos, _neg)
                    proba = output.detach()

                    y_actual.append(target.numpy())
                    y_pred.append(preds.numpy())
                    y_proba.append(proba.numpy())

                    loss = criterion(
                        output, target.float().reshape(target.size()[0], 1, -1)
                    )  # !!! ONLY WORKS FOR CNN !!!

                    if phase == "train":
                        loss.backward()
                        optimizer.step()
                        scheduler.step()  # decay learning rate every few epochs

                running_loss += loss.item() * pep.size(0)
                count += pep.size(0)

            # add epoch labels and predictions to training history
            train_history[phase]["out"]["actual"].append(np.vstack(y_actual))
            train_history[phase]["out"]["pred"].append(np.vstack(y_pred))
            train_history[phase]["out"]["proba"].append(np.vstack(y_proba))

            # calculate metrics for the model at current epoch
            train_history[phase]["metrics"]["loss"].append(running_loss / count)
            train_history[phase]["metrics"] = impepdom.metrics.calculate_metrics(
                train_history
            )[phase]

            epoch_loss = train_history[phase]["metrics"]["loss"][-1]
            epoch_acc = train_history[phase]["metrics"]["acc"][-1]
            epoch_auc = train_history[phase]["metrics"]["auc"][-1]

            wandb.log(
                {
                    f"{phase}_loss": epoch_loss,
                    f"{phase}_acc": epoch_acc,
                    f"{phase}_auc": epoch_auc,
                }
            )

            state_dicts.append(model.state_dict())

            if show_output:
                print(
                    "{} loss: {:.4f} acc: {:.4f} auc: {:.4f}".format(
                        phase, epoch_loss, epoch_acc, epoch_auc
                    )
                )

            # select best model
            if (phase == "val" or not validation) and epoch_auc > best_auc:
                best_auc = epoch_auc
                best_model_wts = copy.deepcopy(model.state_dict())

        if show_output:
            print()  # empty line

    time_elapsed = time.time() - since
    if show_output:
        print(
            "training completed in {:.0f} m {:.4f} s".format(
                time_elapsed // 60, time_elapsed % 60
            )
        )
        print(
            "best {} auc: {:.4f}".format(
                "validation" if validation else "training", best_auc
            )
        )

    model.load_state_dict(best_model_wts)
    return model, train_history, state_dicts
