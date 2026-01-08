from flask import Flask, render_template, jsonify
import Chemistry from
app = Flask(__name__)


#Protein Data
proteins = [
        {
        "id": "VEGF Receptor 2",
        "drugs": [
            {"name": "Sunitinib-", "compound": "placer1"},
            {"name": "Cabozantinib","compound": "placer2"}]
        }
        ,
        {
          "id": "PDGF Receptor beta" ,
          "drugs": [
            {"name": "Sunitinib-", "compound": "placer3"},
            {"name": "Cabozantinib","compound": "placer4"}]
        },
        {
        "id": "mTOR",
        "drugs": [{"name": "Everolimus", "compound": "placer5"}, {"name": "Temsirolimus","compound": "placer1"}]
        }
]
        

@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")

#Getting the protein - Nyosha 
@app.route('/api/proteins', method=['GET'])
def protein_and_drugs_jsonify():
    return jsonify(proteins)
    