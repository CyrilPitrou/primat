# -*- coding: utf-8 -*-
"""
primat.tools
============
Developer-facing generators that derive committed "build product" files from
``primat/config.py``'s ``DEFAULT_PARAMS``/``PARAM_GROUPS``, so keeping them
in sync is "regenerate and commit the diff" instead of hand-editing three
independent copies (see :mod:`primat.tools.gen_param_templates`). Not part of the
public runtime API -- primat itself never imports this package.
"""
