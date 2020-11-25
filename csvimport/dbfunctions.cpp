#include <iostream>
#include <pqxx/pqxx>
#include <string>

#include "dbfunctions.hpp"


void csvimport(pqxx::work& W, std::string datafile, std::string table)
{
  std::string sql_createtbl = "DROP TABLE IF EXISTS " + table + ";"\
      "CREATE TABLE " + table + " (";
    if (table == "dict_estados") {
      sql_createtbl += "id integer, estado text";
    } else if (table == "dict_idades") {
      sql_createtbl += "id integer, idade real";
    } else if (table == "populacao") {
      sql_createtbl += "id integer, populacao text";
    } else if (table == "renda") {
      sql_createtbl += "id integer, fx_1 real, fx_2 real, fx_3 real, fx_4 real,"\
	" fx_5 real, fx_6 real, fx_7 real, share_acum_1 real, share_acum_2 real"\
	", share_acum_3 real, share_acum_4 real, share_acum_5 real, "\
	"share_acum_6 real";
    } else if (table == "idade") {
      sql_createtbl += "id integer, age_0_4 real, age_5_9 real, age_10_14 real, "\
	"age_15_19 real, age_20_24 real, age_25_29 real, age_30_34 real, age_35_"\
	"44 real, age_45_54 real, age_55_64 real, age_65_74 real";
    } else if (table == "lanches_shares" || table == "lanches_precos") {
      sql_createtbl += "j text, I_1 real, I_2 real, I_3 real, I_4 real, I_5 real"\
	", I_6 real, I_7 real, I_8 real, I_9 real, I_10 real, I_11 real, I_12 re"\
	"al, I_13 real, I_14 real, I_15 real, I_16 real, I_17 real, I_18 real, I"\
	"I_1 real, II_2 real, II_3 real, II_4 real, II_5 real, II_6 real, II_7 r"\
	"eal, II_8 real, II_9 real, II_10 real, II_11 real, II_12 real, II_13 re"\
	"al, II_14 real, II_15 real, II_16 real, II_17 real, II_18 real, III_1 r"\
	"eal, III_2 real, III_3 real, III_4 real, III_5 real, III_6 real, III_7 "\
	"real, III_8 real, III_9 real, III_10 real, III_11 real, III_12 real, II"\
	"I_13 real, III_14 real, III_15 real, III_16 real, III_17 real, III_18 r"\
	"eal, IV_1 real, IV_2 real, IV_3 real, IV_4 real, IV_5 real, IV_6 real, "\
	"IV_7 real, IV_8 real, IV_9 real, IV_10 real, IV_11 real, IV_12 real, IV"\
	"_13 real, IV_14 real, IV_15 real, IV_16 real, IV_17 real, IV_18 real, V"\
	"_1 real, V_2 real, V_3 real, V_4 real, V_5 real, V_6 real, V_7 real, V_"\
	"8 real, V_9 real, V_10 real, V_11 real, V_12 real, V_13 real, V_14 real"\
	", V_15 real, V_16 real, V_17 real, V_18 real, VI_1 real, VI_2 real, VI_"\
	"3 real, VI_4 real, VI_5 real, VI_6 real, VI_7 real, VI_8 real, VI_9 rea"\
	"l, VI_10 real, VI_11 real, VI_12 real, VI_13 real, VI_14 real, VI_15 re"\
	"al, VI_16 real, VI_17 real, VI_18 real, VII_1 real, VII_2 real, VII_3 r"\
	"eal, VII_4 real, VII_5 real, VII_6 real, VII_7 real, VII_8 real, VII_9 "\
	"real, VII_10 real, VII_11 real, VII_12 real, VII_13 real, VII_14 real, "\
	"VII_15 real, VII_16 real, VII_17 real, VII_18 real";
    }
    sql_createtbl += ");";

    W.exec(sql_createtbl);

    std::string sql_copy = "COPY " + table + " FROM '" + datafile + "' delimiter ';' NULL as 'N/A' csv header;";
    W.exec(sql_copy);
    std::cout << "Arquivo csv carregado com sucesso: " + table + ".csv" <<\
      std::endl;
}
