#!/usr/bin/env bash

# change cookie consent language
cat _site/site_libs/cookie-consent/cookie-consent.js | perl -npe 's/We use cookies and other tracking technologies to improve your browsing experience on our website, to show you personali[sz]ed content and targeted ads, to analy[sz]e our website traffic, and to understand where our visitors are coming from/We use cookies to analyse our website traffic and to understand where our visitors are coming from/' > tmp.txt

mv tmp.txt _site/site_libs/cookie-consent/cookie-consent.js

# the upload itself
rsync -uvr _site/ ndcn0890@linux.ox.ac.uk:public_html/
