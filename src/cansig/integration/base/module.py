import warnings

from scvi.module.base import BaseModuleClass  # pytype: disable=import-error


_UNTRAINED_WARNING_MESSAGE = "Trying to query inferred values from an untrained model. Please train the model first."


class CanSigBaseModule(BaseModuleClass):
    def __init__(self):
        super(CanSigBaseModule, self).__init__()
        self.is_trained_ = False
        self.history_ = None

    @property
    def is_trained(self) -> bool:
        """Whether the model has been trained."""
        return self.is_trained_

    @is_trained.setter
    def is_trained(self, value):
        self.is_trained_ = value

    @property
    def history(self):
        """Returns computed metrics during training."""
        return self.history_

    def _check_if_trained(self, warn: bool = True, message: str = _UNTRAINED_WARNING_MESSAGE):
        """
        Check if the model is trained.

        If not trained and `warn` is True, raise a warning, else raise a RuntimeError.
        """
        if not self.is_trained_:
            if warn:
                warnings.warn(message)
            else:
                raise RuntimeError(message)
