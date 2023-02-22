#!/usr/bin/env python3
# -*-encoding: utf-8-*-

# created: 20.02.23
# Excusa. Quod scripsi, scripsi.

# by d.zashkonyi

import re
import requests


def get_str_transformations(n, dim=3):
    url = "https://cryst.ehu.es/cgi-bin/cryst/programs/nph-series"
    resp = requests.get(url, params={"gnum": f"{n}"})

    template = '<pre>(.*\n.*\n.*)</pre></b></td>\n<td align=center>(.*)</td></tr>'
    res = re.findall(template, resp.content.decode())
    return res
