# -*- coding: utf-8 -*-
from tests import test_cpd as cpd
from tests import test_mip as mip
from tests import test_mip_reverse as reverse

def main():

    cpd.main()
    mip.main()
    reverse.main()

if __name__ == "__main__":
    main()
