from torch.utils.data import DataLoader


class DataFetch:
    def __init__(self, data):
        self.data = data

    def __getitem__(self, index):
        item = {"tokenized": self.data[index]}
        return item

    def __len__(self):
        return self.data.shape[0]


def molecloader(data, batch_size=64):
    dataset = DataFetch(data)

    return DataLoader(dataset, batch_size=batch_size)
