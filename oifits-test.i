/*
 * TODO:
 *  - check against: <http://jmmc.fr/oival/>
 *  - convert all names to uppercase without trailing spaces
 *  - warn when invalid entries
 */
#include "oifits.i"

func oifits_test
{
  master = oifits_new();
  array_db = oifits_new_array(arrname = "First",
                              frame = "none",
                              arrayx = 0, arrayy = 0, arrayz = 0,
                              tel_name = ["T1", "T2"],
                              sta_name = ["G", "H"],
                              sta_index = [1, 2],
                              diameter = [1.5, 1.5],
                              staxyz = [[ 0.0, 12.0],
                                        [23.0,  0.0],
                                        [ 0.0,  0.0]]);

  target_db = oifits_new_target(target_id = 1,
                                target  = "Binaire",
                                raep0 = 0,
                                decep0 = 0,
                                equinox = 0,
                                ra_err = 0,
                                dec_err = 0,
                                sysvel = 0,
                                veltyp = "OPTICAL",
                                veldef = "UNKNOWN",
                                pmra = 0,
                                pmdec = 0,
                                pmra_err = 0,
                                pmdec_err = 0,
                                parallax = 0,
                                para_err = 0,
                                spectyp = "UNKNOWN");

  wave_db = oifits_new_wavelength(insname = "Lab",
                                  eff_wave = 0.65e-6,
                                  eff_band = 1e-9);

  vis_db = oifits_new_vis(date_obs = "today",
                          arrname = "First",
                          insname = "Lab",
                          target_id = [1,1,1],
                          time = [1,1,1],
                          mjd = [1,1,1],
                          int_time = [1,1,1],
                          visamp = [1,1,1],
                          visamperr = [.1,.1,-.1],
                          visphi = [0,0,0],
                          visphierr = [.1,.1,.1],
                          ucoord = [10,20,10],
                          vcoord = [10,-10,-20],
                          sta_index = [[1,1,1],[2,2,2]]);

  oifits_insert, master, array_db;
  oifits_insert, master, target_db;
  oifits_insert, master, wave_db;
  oifits_insert, master, vis_db;
  oifits_update, master;

  oifits_save, master, "/tmp/test.oifits", overwrite=1;
}

