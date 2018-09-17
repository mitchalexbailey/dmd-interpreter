# _DMD_ Open-access Variant Explorer

### _Contributors: Mitch Bailey, Nicole Miller, Ivo Fokkema_

### Purpose
The purpose of the DMD Open-Access Variant Explorer (DOVE) is to streamline analysis of genetic variants effecting the Dp427m (NM_004006.2) transcript of the _DMD_ gene.

DOVE uses Python and the Django web framework to integrate existing open-access tools to reduce the number od distinct searches needed to analyse a variant, and to expand analysis beyond variants already well described (i.e. novel variants). In addition to predicted consequences of genetic variants, the tool integrates theoretical eligibility for exon skipping or read-through therapies.

### Set up
Clone/download dmd-interpreter then ```cd``` to the dmd-interpreter directory

To activate virtual environment:
```. bin/activate```

To download dependencites (requires python2.7 and corresponding pip):
```pip install -r requirements.txt```

To run locally (suggestions/contributions are welcomed!):
```python manage.py runserver```

Note: For static files to be served locally, ```DEBUG``` will need to be ```TRUE``` in ```settings.py```
