#!/usr/bin/env rbx

require "mysql2"
require "rcsbrest"

puts 1
@con = Mysql2::Client.new(
  host: "localhost",
  username: "root",
  password: "",
  database: "hotspot",
  local_infile: true
)

#method connects and loads to temp folder pdb data
def load_pdb(name)
  BIO
end
