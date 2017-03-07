#!/usr/bin/env slsh
_traceback=1;
set_import_module_path (path_concat (getcwd (),"../src")
			+ ":" + get_import_module_path ());
prepend_to_slang_load_path (path_concat (getcwd (), "../src"));
% putenv ("LD_LIBRARY_PATH=" + path_concat(getcwd(), "../libnovas"));

private define main ()
{
   if (__argc < 2)
     {
	() = fprintf (stderr, "Usage: ./%s <example.sl>\n",
		      path_basename (__argv[0]));
	exit (1);
     }
   variable file = path_concat ("./", __argv[1]);
   __set_argc_argv (__argv[[1:]]);
   () = evalfile (file);
}

main ();
