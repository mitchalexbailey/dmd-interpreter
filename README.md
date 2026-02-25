# _DMD_ Open-access Variant Explorer

### _Contributors: Mitch Bailey, Nicole Miller, Ivo Fokkema_

### Purpose
The purpose of the _DMD_ Open-Access Variant Explorer (DOVE) is to streamline analysis of genetic variants affecting the Dp427m (`NM_004006.2`) transcript of the _DMD_ gene.

DOVE uses Python and the Flask web framework to integrate existing open-access tools to reduce the number of distinct searches needed to analyse a variant, and to expand analysis beyond variants already well described (i.e. novel variants).
In addition to predicted consequences of genetic variants, the tool integrates theoretical eligibility for exon skipping or read-through therapies.

### Set up

- _Requires Python3.x, corresponding pip, virtualenv, Perl_.

- Clone/download dmd-interpreter then `cd` to the dmd-interpreter directory.

- Create and activate the virtual environment:

  ```
  virtualenv DOVE_env --python=python3.11
  . DOVE_env/bin/activate
  ```

- Download all dependencies:

  `pip install -r requirements.txt`

- To run locally (suggestions/contributions are welcomed!):

  `python app.py`

- DOVE can now be tested through your browser. To deactivate:

  ```
  deactivate
  exit
  ```

### Installing DOVE on an Apache2 environment

- _Requires mod-wsgi_.

- Edit Apache's .conf file for the site you wish to add DOVE to.

- Add the following lines to your Apache config (don't forget to change the paths):

  ```
  Alias /static /var/www/dmd-interpreter/interpreter/static
  <Directory /var/www/dmd-interpreter/interpreter/static>
      Require all granted
  </Directory>

  <Directory /var/www/dmd-interpreter>
      <Files interpreter.wsgi>
          Require all granted
      </Files>
  </Directory>

  WSGIDaemonProcess DOVE python-home=/var/www/dmd-interpreter/DOVE_env home=/var/www/dmd-interpreter
  WSGIProcessGroup DOVE
  WSGIScriptAlias /DOVE /var/www/dmd-interpreter/interpreter.wsgi
  ```

- Make sure Apache has read permissions for the application files and static assets.

- Restart Apache to complete the process.

View DOVE installed at the Leiden University Medical Center: http://www.dmd.nl/DOVE
