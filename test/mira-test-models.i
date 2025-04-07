func mira_array_info(str, a)
{
    write, format="%s min = %g, max = %g, avg = %g, rms = %g\n",
        str, min(a), max(a), avg(a), a(*)(rms);
}

func mira_test_models(file, pixelsize=, dims=, img=)
{
    if (is_void(file)) file = get_env("HOME")+"/git/YMiRA/data/data1.oifits";
    if (is_void(pixelsize)) pixelsize = 0.5*MIRA_MILLIARCSECOND;
    if (is_void(img)) img = random((is_void(dims) ? [2,128,128] : dims));
    master = mira_new(file);
    mira_config, master, pixelsize = pixelsize, dims = dimsof(img);
    mira_update, master;
    local x; eq_nocopy, x, master.img_x;
    local y; eq_nocopy, y, master.img_y;
    local u; eq_nocopy, u, master.coords.unique.u;
    local v; eq_nocopy, v, master.coords.unique.v;
    local wave; eq_nocopy, wave, master.coords.unique.wave;
    local band; eq_nocopy, band, master.coords.unique.band;
    write, format="- ufreq: %s\n", (allof(master.coords.ufreq == u/wave) ? "OK" : "ERROR");
    write, format="- vfreq: %s\n", (allof(master.coords.vfreq == v/wave) ? "OK" : "ERROR");
    m = numberof(master.coords.ufreq);
    n1 = numberof(x);
    n2 = numberof(y);

    write, format="- %s:\n", "pixels -> visibilities without bandwidth smearing";
    mira_config, master, xform = "nfft", smearingfunction="none";
    mira_update, master;
    z1 = master.xform(img);
    if (anyof(dimsof(z1) != [2,2,m])) error, "invalid dimensions (NFFT)";
    mira_array_info, "  - Re(V) by NFFT                 ", z1(1,*);
    mira_array_info, "  - Im(V) by NFFT                 ", z1(2,*);
    mira_config, master, xform = "separable", smearingfunction="none";
    mira_update, master;
    z2 = master.xform(img);
    if (anyof(dimsof(z2) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) separable               ", z2(1,*);
    mira_array_info, "  - Im(V) separable               ", z2(2,*);
    mira_array_info, "  - abs. diff.                    ", abs(z2 - z1);
    mira_config, master, xform = "nonseparable", smearingfunction="none";
    mira_update, master;
    z3 = master.xform(img);
    if (anyof(dimsof(z3) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) nonseparable            ", z3(1,*);
    mira_array_info, "  - Im(V) nonseparable            ", z3(2,*);
    mira_array_info, "  - abs. diff.                    ", abs(z3 - z1);
    H = mira_build_separable_pix2vis(x, y, u, v, wave);
    z4 = H(img);
    if (anyof(dimsof(z4) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) alt. model, separable   ", z4(1,*);
    mira_array_info, "  - Im(V) alt. model, separable   ", z4(2,*);
    mira_array_info, "  - abs. diff.                    ", abs(z4 - z1);
    H = mira_build_nonseparable_pix2vis(x, y, u, v, wave);
    z5 = H(img);
    if (anyof(dimsof(z5) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) alt. model, nonseparable", z5(1,*);
    mira_array_info, "  - Im(V) alt. model, nonseparable", z5(2,*);
    mira_array_info, "  - abs. diff.                    ", abs(z5 - z1);

    write, format="- %s:\n", "pixels -> visibilities with bandwidth smearing";
    mira_config, master, xform = "nonseparable", smearingfunction="sinc";
    mira_update, master;
    z7 = master.xform(img);
    if (anyof(dimsof(z7) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) nonseparable            ", z7(1,*);
    mira_array_info, "  - Im(V) nonseparable            ", z7(2,*);
    H = mira_build_nonseparable_pix2vis(x, y, u, v, wave, band, smearingfunction="sinc");
    z8 = H(img);
    if (anyof(dimsof(z8) != [2,2,m])) error, "invalid dimensions (separable)";
    mira_array_info, "  - Re(V) alt. model, nonseparable", z8(1,*);
    mira_array_info, "  - Im(V) alt. model, nonseparable", z8(2,*);
    mira_array_info, "  - abs. diff.                    ", abs(z8 - z7);
}
