
import os

DATA_PATH = "data/tmp/"

DEBUG = True

if os.path.dirname(os.path.realpath(__file__)) == "/home/ivarandknut/checkout_genomic_intervals":
    DATA_PATH = "/var/www/web/data/tmp/"
    #DATA_PATH = "../data/tmp/"

