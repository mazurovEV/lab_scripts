
# coding: utf-8

# In[ ]:


import flask
import json
import pandas as pd
import numpy as np
from neo4j import GraphDatabase
from pathlib import Path
from flask_cors import CORS
import json
from datetime import datetime

app = flask.Flask(__name__, static_folder='dist', static_url_path="/dist")
#app = flask.Flask(__name__)
CORS(app)

# app.config["DEBUG"] = True

driver = GraphDatabase.driver("bolt://localhost:7687")

def get_db():
    if not hasattr(flask.g, 'neo4j_db'):
        flask.g.neo4j_db = driver.session()
    return flask.g.neo4j_db


@app.teardown_appcontext
def close_db(error):
    if hasattr(flask.g, 'neo4j_db'):
        flask.g.neo4j_db.close()

@app.route("/")
@app.route("/browser")
@app.route("/search")
def catch_all():
    return app.send_static_file("index.html")


# In[ ]:


@app.route('/api/v1/dashboard', methods=['GET'])
def browser():
    response = {
        "table": {
            "all_counts": 11,
            "data": [
                {"HM": "H3K27ac",
                "number of lncRNAs with corr": 2031,
                "number of peaks with corr": 841674,
                "number of genes associated with peaks": 50788,
                "number of tissues": 49},
                
                {"HM": "H3K27me3",
                "number of lncRNAs with corr": 1401,
                "number of peaks with corr": 198210,
                "number of genes associated with peaks": 43544,
                "number of tissues": 50},
                
                {"HM": "H3K36me3",
                "number of lncRNAs with corr": 1321,
                "number of peaks with corr": 342718,
                "number of genes associated with peaks": 49875,
                "number of tissues": 52},
                
                {"HM": "H3K4me1",
                "number of lncRNAs with corr": 1954,
                "number of peaks with corr": 906110,
                "number of genes associated with peaks": 53258,
                "number of tissues": 51},
                
                {"HM": "H3K4me2",
                "number of lncRNAs with corr": 1031,
                "number of peaks with corr": 447751,
                "number of genes associated with peaks": 41151,
                "number of tissues": 19},
                
                {"HM": "H3K4me3",
                "number of lncRNAs with corr": 1752,
                "number of peaks with corr": 738667,
                "number of genes associated with peaks": 46337,
                "number of tissues": 59},
                
                {"HM": "H3K9ac",
                "number of lncRNAs with corr": 1104,
                "number of peaks with corr": 78177,
                "number of genes associated with peaks": 23634,
                "number of tissues": 19},
                
                {"HM": "H3K9me3",
                "number of lncRNAs with corr": 1079,
                "number of peaks with corr": 250631,
                "number of genes associated with peaks": 38708,
                "number of tissues": 50},
                
                {"HM": "H3K79me2",
                "number of lncRNAs with corr": 802,
                "number of peaks with corr": 102280,
                "number of genes associated with peaks": 32919,
                "number of tissues": 20},
                
                {"HM": "H4K20me1",
                "number of lncRNAs with corr": 841,
                "number of peaks with corr": 61974,
                "number of genes associated with peaks": 22777,
                "number of tissues": 19}
            ]
        }
    }
    
    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response


# In[ ]:

def _search(tx, query):
    result = tx.run(query)
    return [{"Histone Modification": record["hm"][2], "lncRNA": record["lnc"], "Peak Id": record["name"], "Chr": record["chrom"], "Start": str(record["start"]), "End": str(record["end"]), "Gene": str(record["gene"]), "Corr": str(record["corr"])} for record in result]

def _manhattan(tx, query):
    result = tx.run(query)
    t = [{"Histone Modification": record["hm"][2], "Chr": record["chr"], "start": int(record["start"]), "end": int(record["end"]), "corr": str(record["corr"])} for record in result]
    
    df = pd.DataFrame(t)
    df['center'] = df["start"] + np.ceil((df["end"] - df["start"])//2)
    df.to_csv("/home/mazurovev/site_log/before_mn.tsv", sep="\t", index=None)
    df = df.apply(lambda row: [row["Histone Modification"], row['Chr'], row["center"]/(10**(len(str(row["center"])) - 3)), row['corr']], axis=1, result_type='expand')
    df.to_csv("/home/mazurovev/site_log/mn.tsv", sep="\t", index=None)
    res = {}
    for hm in df[0].unique():
        tmp = df[df[0] == hm]
        
        res[hm] = []
        for _ in range(24):
            res[hm].append([]) 
        for i, row in tmp.iterrows():
            c = row[1].split("chr")[1]
            if "chr" in row[1] and (c >= "1" or c <= "22" or c =="X" or c == "Y"):
                if c == "X":
                    number = 23
                elif c == "Y":
                    number = 24
                else:
                    number = int(c)
                
                res[hm][number - 1].append({"x": (number - 1)*10 + row[2], "y": row[3]})
    with open(Path("/home/mazurovev/site_log/chr.txt"), "w") as f:
        f.write("\n\n" + str(res) + "\n")    
    return res

def _search_other_peaks(tx, query):
    result = tx.run(query)
    return [{"Histone Modification": record["hm"][2], "Chr": record["chrom"], "Start": record["start"], "End": record["end"], "Peak Id": record["name"]} for record in result]

def _expression(tx, query):
    result = tx.run(query)
    return [{"label": record["tissue"], "value": str(round(float(record["expression"]), 2))} for record in result]

def _corr(tx, query):
    result = tx.run(query)
    return [{"signal": str(round(float(record["signal"]), 2)), "expression": str(round(float(record["expression"]), 2)), "tissue": record["tissue"]} for record in result]

def _count(tx, query):
    result = tx.run(query)
    return [record["count"] for record in result][0]

def sort_chroms(chrom):
    return int(chrom.split('chr')[1] if chrom.split('chr')[1]!='X' and chrom.split('chr')[1] != 'Y' else 40 if chrom.split('chr')[1]=='X' else 41)

def _count_chroms(tx, query):
    result = tx.run(query)
    tmp = [{"label": record["p.chrom"], "value": record["count(p)"]} for record in result]
    return sorted(tmp, key=lambda x: sort_chroms(x['label']))

def _get_cypher_query(json_data, limit=True):
    hm_queries = []
    count_queries = []
    result = "Return labels(p) AS hm, lnc.name AS lnc, p.chrom AS chrom, p.start AS start, p.end AS end, gene.name AS gene, c.corr AS corr, p.name AS name"
    result_count = "Return p"
    for hm in json_data["hm"]:
        base_query = "Match (p: " + hm + ")-[c:CORR]-(lnc:lncRNA)"
        lnc_where = []
        coords_where = ""
        corrs_where = ""
        gene_match = ""
        skip_limit = ""
        
        if json_data["lncrna"]:
            for l in json_data["lncrna"]:
                if "ENSG" in l:
                    lnc_where.append("lnc.ensembl_id='" + l + "'")
                else:
                    lnc_where.append("lnc.name='" + l + "'")
        
            lnc_where = "(" + " OR ".join(lnc_where) + ")"
        
        coords = []
        if json_data["coords"]:
            for i, c in enumerate(json_data["coords"]):
                coords.append("p.chrom=" + "'" + c[0] + "'" + " AND p.start >= " + c[1] + " AND p.end <= " + c[2])
        coords_where ="(" + " OR ".join(coords) + ")" if coords else ""
        
        corrs = []
        if json_data["thresholds_choisen"]:
            if json_data["thresholds_choisen"][0]:  # пороги на + корреляцию
                corrs.append("c.corr >= " + json_data["tresholds"][0][0] + " AND c.corr <= " + json_data["tresholds"][0][1])
            if json_data["thresholds_choisen"][1]:  # пороги на - корреляцию
                corrs.append("c.corr >= " + json_data["tresholds"][1][0] + " AND c.corr <= " + json_data["tresholds"][1][1])
        corrs_where = "(" + " OR ".join(corrs) + ")" if corrs else ""
        
        if json_data["genes"]:
            gene_match = "Match (p)-[:IN]-(gene:Gene)"
            gene_where = []
            for l in json_data["genes"]:
                if "ENSG" in l:
                    gene_where.append("gene.ensembl_id='" + l + "'")
                else:
                    gene_where.append("gene.name='" + l + "'")
        
            gene_where = "(" + " OR ".join(gene_where) + ")"
            
            gene_match = gene_match + " WHERE " + gene_where
        else:
            gene_match = "Optional Match (p)-[:IN]-(gene:Gene)"
        
        where_clause = []
        if lnc_where:
            where_clause.append(lnc_where)
        if corrs_where:
            where_clause.append(corrs_where)
        if coords_where:
            where_clause.append(coords_where)
        
        hm_queries.append(base_query + (" WHERE " if where_clause else "") + " AND ".join(where_clause) + " " + gene_match + " " + result + " " + skip_limit)
        
        count_queries.append(base_query + (" WHERE " if where_clause else "") + " AND ".join(where_clause) + " " + gene_match + " " + result_count)
        
    page = int(json_data["page"])
    page_count = int(json_data["page_count"])
    skip_limit = "SKIP " + str((page - 1) * page_count) + " LIMIT " + str(page_count)
    
    return "CALL { " + " UNION ".join(hm_queries) + " } " + "Return hm, lnc, chrom, start, end, gene, corr, name" + (" " + skip_limit if limit else ""), count_queries

@app.route('/api/v1/search/results', methods=['GET', 'POST'])
def search():
    json_data = flask.request.get_json(force=True)
    
    db = get_db()
    
    cypher, count_queries = _get_cypher_query(json_data)
   
    cypher_count = "CALL { " + " UNION ".join(count_queries) + " } " + "Return count(p) AS count"
    
    result = db.read_transaction(_search, cypher)
    # count_result = db.read_transaction(_count, cypher_count)
    if Path("/home/mazurovev/site_log/ZBTB33.tsv").is_file():
        pd.DataFrame(result).to_csv("/home/mazurovev/site_log/ZBTB33.tsv", index=None, mode="a", sep="\t", header="False")
    else:
        pd.DataFrame(result).to_csv("/home/mazurovev/site_log/ZBTB33.tsv", index=None, sep="\t")
    response = {
        "title": "Search result",
        "table": {
            #"all_counts": count_result,
            "data": result
        }
    }
    
    # with open(Path("/home/mazurovev/site_log/log.txt"), "a") as f:
    #    f.write("\n\n" + str(datetime.now()) + "\n")
    #    f.write(json.dumps(json_data))
    #    f.write("\n")
    #    f.write(cypher)
    #    for i in result:
    #        f.write(str(i) + "\n")
    #    f.write(json.dumps(response))
    #    f.write("\n\n")
    #    f.write(cypher_count)
        # f.write("count: " + str(count_result))

    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    response.headers.add("Access-Control-Allow-Headers", "*")

    return response


# In[26]:


@app.route('/api/v1/info/modification', methods=['GET'])
def modification():
    modification = flask.request.args.get('hm')
    page = int(flask.request.args.get('page'))
    page_count = int(flask.request.args.get('page_count'))
    
    cypher = 'match (p:' + modification + ') where p.chrom in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"] return p.chrom, count(p);'
    
    table_cypher = 'match (p:' + modification + ')-[c:CORR]-(lnc:lncRNA) Optional match (p)-[:IN]-(g:Gene) Return labels(p) AS hm, lnc.name AS lnc, p.chrom AS chrom, p.start AS start, p.end AS end, g.name AS gene, c.corr AS corr, p.name AS name SKIP ' + str(page*page_count) + ' LIMIT ' + str(page_count)
    
    cypher_count = 'match (p:' + modification + ')-[c:CORR]-(lnc:lncRNA) Optional match (p)-[:IN]-(g:Gene) Return count(p) AS count;'
    
    db = get_db()
    
    response = {
        "chart": {
            "title": modification + " peaks distribution",
            "elements": db.read_transaction(_count_chroms, cypher)
        },
        "table": {
            # "all_counts": db.read_transaction(_count, cypher_count),
            "title": "peaks of " + modification + " modification",
            "data": db.read_transaction(_search, table_cypher)
        }
     }

    
    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response


# In[29]:


@app.route('/api/v1/info/lncrna', methods=['GET'])
def lncrna():
    lncrna = flask.request.args.get('lncrna')

    page = int(flask.request.args.get('page'))
    page_count = int(flask.request.args.get('page_count'))
    
    cypher = 'match (lnc:lncRNA {name: "' + lncrna + '"})-[e:EXPRESS_IN]-(t:Tissue) where e.expression<>0 return distinct t.name AS tissue, avg(e.expression) AS expression order by tissue'
    
    table_cypher = 'match (lnc:lncRNA {name: "' + lncrna + '"})-[c:CORR]-(p:Peak) Optional match (p)-[:IN]-(g:Gene) Return labels(p) AS hm, lnc.name AS lnc, p.chrom AS chrom, p.start AS start, p.end AS end, g.name AS gene, c.corr AS corr, p.name AS name SKIP ' + str(page*page_count) + ' LIMIT ' + str(page_count)
    
    cypher_count = 'match (lnc:lncRNA {name: "' + lncrna + '"})-[c:CORR]-(p:Peak) Optional match (p)-[:IN]-(g:Gene) Return count(p) AS count;'
    
    cypher_manhettan = 'match (p:Peak)-[c:CORR]-(lnc: lncRNA {name: "' + lncrna + '"}) where c.corr > 0.6 or c.corr < -0.6 return labels(p) AS hm, p.chrom AS chr, p.start AS start, p.end AS end, c.corr AS corr'
    
    db = get_db()
    
    response = {
        "chart": {
            "title": lncrna + " expression",
            "elements": db.read_transaction(_expression, cypher)
        },
        "barplot": {
            "labels": ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", "Chr22", "ChrX", "ChrY"],
            "elements": db.read_transaction(_manhattan, cypher_manhettan)
        },
        "table": {
            # "all_counts": db.read_transaction(_count, cypher_count),
            "title": "peaks of " + lncrna,
            "data": db.read_transaction(_search, table_cypher)
        }
     }
    
    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response


# In[31]:
@app.route('/api/v1/download', methods=['POST', 'GET'])
def download():
    #json_data = json.loads(flask.request.args.get('json'))
    with open(Path("/home/mazurovev/site_log/test_dw.txt"), "a") as f:
        f.write("Работает!")
    json_data = flask.request.get_json(force=True)
        
    db = get_db()
    
    cypher, _ = _get_cypher_query(json_data, limit=False)
    result = db.read_transaction(_search, cypher)
    
    def generate():
        res = [{"Header_1": "Histone Modification"}, {"Header_2": "lncRNA"}, {"Header_3": "Chr"}, {"Header_4": "Start"}, {"Header_5": "End"}, {"Header_6": "gene_id"}, {"Header_7": "gene_name"}, {"Header_8": "correlation"}, {"Header_9": "peak_id"}] + result
        for d in res:
            yield bytes('\t'.join(d.values()) + '\n', encoding="utf-8")
    
    pd.DataFrame(result).to_csv("/home/mazurovev/site_log/ZBTB33.tsv", index=None, sep="\t")
    
    resp = flask.Response(generate(), status=200, mimetype='text/csv', direct_passthrough=True)
    
    resp.headers["Content-Disposition"] = "attachment; filename=export.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp

@app.route('/api/v1/info/gene', methods=['GET'])
def gene():
    gene = flask.request.args.get('gene')
    
    page = int(flask.request.args.get('page'))
    page_count = int(flask.request.args.get('page_count'))
    
    other_page = int(flask.request.args.get('other_page'))
    other_page_count = int(flask.request.args.get('other_page_count'))
    
    cypher = 'match (g:Gene {name: "' + gene + '"})-[e:EXPRESS_IN]-(t:Tissue) where e.expression<>0 return distinct t.name AS tissue, avg(e.expression) AS expression order by tissue'
    
    table_cypher = 'match (g: Gene {name: "' + gene + '"})-[:IN]-(p:Peak)-[c:CORR]-(lnc:lncRNA) Return labels(p) AS hm, lnc.name AS lnc, p.chrom AS chrom, p.start AS start, p.end AS end, g.name AS gene, c.corr AS corr, p.name AS name SKIP ' + str(page*page_count) + ' LIMIT ' + str(page_count)
    
    other_peaks_cypher = 'match (g: Gene {name: "' + gene + '"})-[:IN]-(p:Peak) return labels(p) AS hm, p.name AS name, p.chrom AS chrom, p.start AS start, p.end AS end order by p.start, p.end SKIP ' + str(other_page*other_page_count) + ' LIMIT ' + str(other_page_count)
    
    cypher_count = 'match (g: Gene {name: "' + gene + '"})-[:IN]-(p:Peak)-[c:CORR]-(lnc:lncRNA) Return count(p) AS count;'
        
    cypher_other_count = 'match (g: Gene {name: "' + gene + '"})-[:IN]-(p:Peak) return count(p) AS count;'
    
    db = get_db()
    
    response = {
        "chart": {
            "title": gene + " expression",
            "elements": db.read_transaction(_expression, cypher)
        },
        "table": {
            "title": "peaks with correlations",
            #"all_counts": db.read_transaction(_count, cypher_count),
            "data": db.read_transaction(_search, table_cypher)
        },
        "other_peaks_table": {
            "title": "other peaks annotated by this gene",
            #"all_counts": db.read_transaction(_count, cypher_other_count),
            "data": db.read_transaction(_search_other_peaks, other_peaks_cypher)
        }
     }
    
    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response


# In[34]:


@app.route('/api/v1/info/corr', methods=['GET'])
def corr():
    lncrna = flask.request.args.get('lncrna')
    hm = flask.request.args.get('hm')
    peak_id = flask.request.args.get('peak_id')
    
    cypher = 'match (p:' + hm + ' {name: "' + peak_id + '"})-[c:CHIPSEQ_IN]-(t:Tissue)-[e:EXPRESS_IN]-(lnc:lncRNA {name: "' + lncrna + '"}) return t.name AS tissue, c.signal AS signal, avg(e.expression) AS expression'
    
    db = get_db()
    res = db.read_transaction(_corr, cypher)
    
    res_chart = [{"x": r["signal"], "y": r["expression"]} for r in res]
    
    res_signal  = {r["tissue"] : r["signal"] for r in res}
    res_expression  = {r["tissue"] : r["expression"] for r in res}
    res_table = [res_signal, res_expression]
    
    response = {"response": {
                    "title": hm + " peak " + peak_id + " lncRNA " + lncrna + " correlation plot",
                    "chart": {
                        "points": res_chart
                    },
                    "table": {
                        "title": "Correlation data",
                        "data": res_table
                    }
    }}
    
    response = flask.Response(json.dumps(response, ensure_ascii=False), status=200, mimetype='application/json')
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response

