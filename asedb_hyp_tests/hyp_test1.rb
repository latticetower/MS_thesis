#!/usr/bin/env rbx

require "mysql2"

puts 1
@con = Mysql2::Client.new(
  host: "localhost",
  username: "root",
  password: "",
  database: "hotspot",
  local_infile: true
)
