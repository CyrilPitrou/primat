# `primat.network_data` — networks and rate tables

Loading a network: which reactions it holds, the rate table behind each, the
master T9 grid they are all resampled onto, and the A/Z conservation checks
that reject a malformed addition at load time rather than mis-integrating it.

`UpdateNuclearRates` documents the `custom_network=` schema
(`{"removed": [...], "replaced": {...}, "added": {...}}`) that the GUI's
exported archives and the API share — see {doc}`../howto/custom-networks`.

```{eval-rst}
.. automodule:: primat.network_data
   :members:
   :undoc-members:
   :show-inheritance:
```
