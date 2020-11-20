#include <iostream>
#include <pqxx/pqxx>
#include <string>

#include "dbfunctions.hpp"

using std::string;

void csvimport_anac(pqxx::work& W, string datadir, string year, string month)
{
    string table = "csv_" + year + month;
    string sql_createtbl = "DROP TABLE IF EXISTS " + table + ";"
        "CREATE TABLE " + table + " ("
        "ano         integer,"
        "mes         integer,"
        "empresa     text,"
        "origem      text,"
        "destino     text,"
        "tarifa      real,"
        "assentos    integer"
        ");";

    W.exec(sql_createtbl);

    string sql_copy = "COPY " + table + " FROM '" + datadir + year + month + ".csv' DELIMITER ';' CSV HEADER;";
    W.exec(sql_copy);
    string drop_columns = "ALTER TABLE " + table + " DROP COLUMN ano, DROP COLUMN mes;";
    W.exec(drop_columns);
    std::cout << "Arquivo csv carregado com sucesso: anac " << year << month << std::endl;
}

void csvimport_mktsinfo(pqxx::work& W, std::string datafile)
{
    string table = "mktsinfo";
    string sql_createtbl = "DROP TABLE IF EXISTS " + table + ";"
        "CREATE TABLE " + table + " ("
        "id                   integer,"
        "mercado              text,"
        "municipio_origem     text,"
        "municipio_destino    text,"
        "dist                 real,"
        "pop_origem           real,"
        "pop_destino          real,"
        "media_g_pop          real,"
        "media_pop            real"
        ");";

    W.exec(sql_createtbl);

    string sql_copy = "COPY " + table + " FROM '" + datafile + "' delimiter ',' csv header;";
    W.exec(sql_copy);
    std::cout << "Arquivo csv carregado com sucesso: merc_stats" << std::endl;
}

void csvimport_ipca(pqxx::work& W, std::string datafile)
{
    string table = "ipca";
    string sql_createtbl = "DROP TABLE IF EXISTS " + table + ";"
        "CREATE TABLE " + table + " ("
        "date                 text,"
        "index                real"
        ");";

    W.exec(sql_createtbl);

    string sql_copy = "COPY " + table + " FROM '" + datafile + "' delimiter ';';";
    W.exec(sql_copy);
    std::cout << "Arquivo csv carregado com sucesso: ipca" << std::endl;
}

void csvimport_instruments(pqxx::work& W, std::string datafile1)
{
    string table = "instruments";
    string sql_createtbl = "DROP TABLE IF EXISTS " + table + ";"
        "CREATE TABLE " + table + " ("
        "date                 text,"
        "querosene            real"
        ");";

    W.exec(sql_createtbl);

    string sql_copy1 = "COPY " + table + " FROM '" + datafile1 + "' delimiter ';';";
    W.exec(sql_copy1);
    std::cout << "Tabela carregada com sucesso: instruments" << std::endl;
}
