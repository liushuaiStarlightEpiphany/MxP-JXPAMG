
#get svn reversion
function get_svn_reversion()
{
  src_dir  = $1
  svn_info = $(LANG=C svn info $src_dir)
  if [ $? != 0 ]; then
    return 0
  fi

  reversion=$(LANG=C svn info ${src_dir} | grep "Revision" | cut -d " " -f 2)
  return ${reversion}
}
export -f get_svn_reversion
# e.g: get_svn_reversion <svn_src_directory>


