#include <iostream>
#include <pqxx/pqxx>
#include <string>
#include <vector>

#include "dbfunctions.hpp"

using std::string;
using std::vector;


int main(int argc, char* argv[])
{
    try {
        pqxx::connection C("dbname = aviacao user = postgres password = passwd"\
               " hostaddr = 127.0.0.1 port = 5432");
        if (C.is_open()) {
            std::cout << "Opened database successfully: " << C.dbname() << \
                    std::endl;
        } else {
            std::cout << "Can't open database" << std::endl;
            return 1;
        }

        // csv import from ANAC
        if (argc > 1 && std::strcmp(argv[1], "anac") == 0) {
            pqxx::work W0(C);
            vector<string> anos = {"2002", "2003", "2004", "2005", "2006", \
                    "2007", "2008", "2009", "2010", "2011", "2012", "2013", \
                    "2014", "2015", "2016", "2017", "2018", "2019"};
            vector<string> meses = {"01", "02", "03", "04", "05", "06", "07", \
                    "08", "09", "10", "11", "12"};
            string datadir = "/var/lib/postgresql/data/pgdata/databases/anac/";
            for (string ano : anos) {
                for (string mes : meses) {
                    /* 2020 csvs do not look correct to this moment, skipping
                    if (ano == "2020" && !(mes == "01" || mes == "02"))
                        continue;
                    */
                    csvimport_anac(W0, datadir, ano, mes);
                }
            }
            W0.commit();
        }

        // csv import mktsinfo
        if (argc > 1 && std::strcmp(argv[1], "mktsinfo") == 0) {
            pqxx::work W1(C);
            string datafile = "/var/lib/postgresql/data/pgdata/databases/mktsinfo/merc_stats.csv";
            csvimport_mktsinfo(W1, datafile);
            W1.commit();
        }

        // csv import IPCA
        if (argc > 1 && std::strcmp(argv[1], "ipca") == 0) {
            pqxx::work W2(C);
            string datafile = "/var/lib/postgresql/data/pgdata/databases/ipca/ipca.csv";
            csvimport_ipca(W2, datafile);
            W2.commit();
        }

        // csv import instruments
        if (argc > 1 && std::strcmp(argv[1], "instruments") == 0) {
            pqxx::work W3(C);
            string datafile1 = "/var/lib/postgresql/data/pgdata/databases/instruments/indiceQueroseneAviacao.csv";
            csvimport_instruments(W3, datafile1);
            W3.commit();
        }

        std::cout << "Database operation successfull" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
