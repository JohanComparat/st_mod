# !/usr/bin/env python
import sys, os
import logging
from astropy.table import Table

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    handlers=[
        logging.FileHandler('move_eSASS.log'),
        logging.StreamHandler(sys.stdout),
    ],
)
os.environ['UCHUU']='/home/idies/workspace/erosim/Uchuu'
os.environ['GIT_STMOD']='/home/idies/workspace/erosim/software/st_mod'
os.environ['GIT_STMOD_DATA']='/home/idies/workspace/erosim/software/st_mod_data'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

eRASSn = 's4'
for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    field_id = str(sky_tile['SRVMAP']).zfill(6)
    indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", field_id, eRASSn + '_eSASS')
    outdir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", field_id, eRASSn + '_eSASS_pre_2026_06_01')
    cmd = 'mv ' + indir + ' ' + outdir
    logging.info('CMD: %s', cmd)
    rc = os.system(cmd)
    if rc != 0:
        logging.warning('FAILED (exit %d): %s', rc, cmd)
    else:
        logging.info('OK: %s', cmd)

eRASSn = 's5'
for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    field_id = str(sky_tile['SRVMAP']).zfill(6)
    indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", field_id, eRASSn + '_eSASS')
    outdir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", field_id, eRASSn + '_eSASS_pre_2026_06_01')
    cmd = 'mv ' + indir + ' ' + outdir
    logging.info('CMD: %s', cmd)
    rc = os.system(cmd)
    if rc != 0:
        logging.warning('FAILED (exit %d): %s', rc, cmd)
    else:
        logging.info('OK: %s', cmd)
