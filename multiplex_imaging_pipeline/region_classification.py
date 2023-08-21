import numpy as np
import pytorch_lightning as pl
import torchvision.transforms.functional as TF
import torch
import torch.nn as nn
from torchvision.transforms import Compose, ColorJitter, Normalize, RandomAffine, RandomResizedCrop
from torch.utils.data import Dataset, DataLoader
from einops import rearrange
from sklearn.metrics import accuracy_score, f1_score
from timm import create_model
from pytorch_lightning.callbacks import ModelCheckpoint






class RegionImgTransform(object):
    def __init__(self, p=.8, brightness=.2, contrast=.2, saturation=.2, hue=.2,
                 means=None, stds=None):
        
        self.affine_transforms = Compose([
            RandomAffine(180),
            RandomResizedCrop((256, 256), scale=(.6, 1.)),
        ])
        
        self.color_transforms = Compose([
            ColorJitter(brightness=brightness, contrast=contrast, saturation=saturation, hue=hue),
        ])
        
        if means is not None:
            self.normalize = Normalize(means, stds)
        else:
            self.normalize = nn.Identity()
 
        self.p = p
    
    def __call__(self, he, mask):
        """
        """
        if torch.rand(size=(1,)) < self.p:
            he = self.color_transforms(he)
            
            combined = torch.concat((he, mask))
            combined = self.affine_transforms(combined)
            
            he, mask = combined[:-1], combined[-1:]
            
        he = self.normalize(he)
        
        return he, mask


class RegionDataset(Dataset):
    def __init__(self, region_to_cls_imgs, transform=None):
        self.region_to_cls_imgs = region_to_cls_imgs
        self.transform = transform
        self.keys = sorted(self.region_to_cls_imgs.keys())
        
        self.labels = sorted(set([d['label'] for reg, d in region_to_cls_imgs.items()]))
        self.y = [self.labels.index(self.region_to_cls_imgs[k]['label']) for k in self.keys]
        self.y = torch.nn.functional.one_hot(torch.tensor(
            self.y)).to(torch.float32)

    def __len__(self):
        return self.y.shape[0]

    def __getitem__(self, idx):
        region_id = self.keys[idx]

        rgb = self.region_to_cls_imgs[region_id]['rgb']
        mask = self.region_to_cls_imgs[region_id]['mask']
        mask = mask.to(torch.float32)

        if self.transform:
            rgb, mask = self.transform(rgb, mask)

        return {
            'mask': mask,
            'rgb': rgb,
            'y': self.y[idx],
        }
    

class RegionClassifier(nn.Module):
    def __init__(
            self,
            n_classes,
            n_features=100,
            n_channels=3,
            backbone='resnet34',
    ):
        super().__init__()
        self.rgb_encoder = create_model(backbone, in_chans=n_channels)
        self.mask_encoder = create_model(backbone, in_chans=1)
        
        self.encoded_downsample = nn.Linear(1000, n_features)
        
        self.head = nn.Sequential(
            nn.Linear(n_features, n_features),
            nn.BatchNorm1d(n_features),
            nn.ReLU(),
            nn.Dropout(.2),
            nn.Linear(n_features, n_features),
            nn.BatchNorm1d(n_features),
            nn.ReLU(),
            nn.Dropout(.2),
            nn.Linear(n_features, n_classes),
            nn.Softmax(dim=-1)
        )
        
        self.loss = nn.CrossEntropyLoss()   
        
    def calculate_loss(self, y_pred, y_true):
        return self.loss(y_pred, y_true)

    def forward(self, rgb, mask):
        rgb_encoded = self.rgb_encoder(rgb)
        mask_encoded = self.mask_encoder(mask)
        encoded = rgb_encoded + mask_encoded
        encoded = self.encoded_downsample(encoded)

        probs = self.head(encoded)
        
        return probs


class ModelLightning(pl.LightningModule):
    def __init__(self, model, lr=1e-4):
        super().__init__()
        
        self.model = model
        self.lr = lr
        
#         self.save_hyperparameters(ignore=['model'])
        self.save_hyperparameters()

    def training_step(self, batch, batch_idx):
        rgb, mask, y = batch['rgb'], batch['mask'], batch['y']
        probs = self.model(rgb, mask)
        loss = self.model.calculate_loss(probs, y)
        acc = accuracy_score(y.argmax(dim=-1).clone().detach().cpu().numpy(),
                             probs.argmax(dim=-1).clone().detach().cpu().numpy())
        f1 = f1_score(y.argmax(dim=-1).clone().detach().cpu().numpy(),
                      probs.argmax(dim=-1).clone().detach().cpu().numpy(), average='weighted')
        
        self.log_dict({'train/loss': loss, 'train/acc': acc, 'train/f1': f1},
                      on_step=False, on_epoch=True, prog_bar=True)
        
        
        result = {'loss': loss, 'probs': probs, 'y': y}
        return result
    
    def validation_step(self, batch, batch_idx):
        rgb, mask, y = batch['rgb'], batch['mask'], batch['y']
        probs = self.model(rgb, mask)
        loss = self.model.calculate_loss(probs, y)
        acc = accuracy_score(y.argmax(dim=-1).clone().detach().cpu().numpy(),
                             probs.argmax(dim=-1).clone().detach().cpu().numpy())
        f1 = f1_score(y.argmax(dim=-1).clone().detach().cpu().numpy(),
                      probs.argmax(dim=-1).clone().detach().cpu().numpy(), average='weighted')
        
        self.log_dict({'val/loss': loss, 'val/acc': acc, 'val/f1': f1},
                      on_step=False, on_epoch=True, prog_bar=True)
        
        
        result = {'loss': loss, 'probs': probs, 'y': y}
        return result
    
    def prediction_step(self, batch, batch_idx):
        rgb, mask, y = batch['rgb'], batch['mask'], batch['y']
        probs = self.model(rgb, mask)

        result = {'y_pred': probs.argmax(dim=-1), 'probs': probs, 'y_true': y}
        return result

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        return optimizer
    
    def forward(self, batch):
        rgb, mask, y = batch['rgb'], batch['mask'], batch['y']
        return self.model(rgb, mask)
    
