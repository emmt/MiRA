func mira_plugin_example_init(nil) {
  write, "PLUGIN EXAMPLE";
  return mira_new_plugin(options=_lst("\nExample plugin options:",
                                      _lst("example_option", 0, "NUMBER", OPT_INTEGER,
                                           "Some unused value")),
                         parse_options = example_parse_options,
                         tweak_visibilities = example_tweak_visibilities);
}

func example_parse_options(plugin, opt)
{
  h_set, plugin, param = opt.example_option, warn=1n;
}

func example_tweak_visibilities(master, vis)
{
  plugin = mira_plugin(master);
  if (plugin.warn) {
    inform, "EXAMPLE: param = %d", plugin.param;
    h_set, plugin, warn=0n;
  }
  return vis;
}
