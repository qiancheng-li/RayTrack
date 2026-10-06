import os
import tempfile
from pathlib import Path
import shapely
import osmnx as ox
from utils import *
import open3d as o3d
from tqdm import tqdm
from pyproj import Transformer
import matplotlib.pyplot as plt
from triangle import triangulate
from PIL import Image, ImageDraw
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
from shapely.geometry import shape, Polygon

import logging

osm_server_addr = "https://overpass-api.de/api/interpreter"
points = [(-122.312299126737, 47.65074792715546), (-122.312299126737, 47.65690207846001), (-122.30187970485449, 47.65690207846001), (-122.30187970485449, 47.65074792715546)]
ground_scale = 1.5


# 1. setup server and transform

#osm server
ox.settings.overpass_url = osm_server_addr #osm api
ox.settings.overpass_rate_limit = False #false for personal use

#gps projection
EPSG_code = get_utm_epsg_code_from_gps(points[0][0], points[0][1]) #latlon to xy
print(EPSG_code)
to_utm = Transformer.from_crs("EPSG:4326", EPSG_code, always_xy=True) #4326 for gps
to_gps = Transformer.from_crs(EPSG_code, "EPSG:4326", always_xy=True)



# 2. material constant

# Output paths
project_root = Path(__file__).resolve().parent.parent
asset_dir = Path(os.environ.get(
    "RAYLOC_SCENE_OUTPUT_DIR",
    project_root / "assets" / "raytrack_scene",
))
mesh_dir = asset_dir / "test_mesh"
mesh_dir.mkdir(parents=True, exist_ok=True)

#rendering defauult
spp_default = 4096 #samples per pixel
resx_default = 1024
resy_default = 768

camera_settings = {
    "rotation": (0, 0, -90),  # look down
    "fov": 42.854885 #field of view
}

material_colors = {
    "mat-itu_concrete": (0.539479, 0.539479, 0.539480),
    "mat-itu_marble": (0.701101, 0.644479, 0.485150),
    "mat-itu_metal": (0.219526, 0.219526, 0.254152),
    "mat-itu_wood": (0.043, 0.58, 0.184),
    "mat-itu_medium_dry_ground": (0.91, 0.569, 0.055),
}



# 3. xml header

#default
scene = ET.Element("scene", version = "2.1.0")
ET.SubElement(scene, "default", name = "spp", value = str(spp_default))
ET.SubElement(scene, "default", name = "resx", value = str(resx_default))
ET.SubElement(scene, "default", name = "resy", value = str(resy_default))

#integrator
integrator = ET.SubElement(scene, "integrator", type="path")
ET.SubElement(integrator, "integer", name="max_depth", value="12")

#material
for material, value in material_colors.items():
    bsdf = ET.SubElement(scene, "bsdf", type = "twosided", id = material)
    bsdf_diffuse = ET.SubElement(bsdf, "bsdf", type = "diffuse")
    ET.SubElement(bsdf_diffuse, "rgb", value = f"{value[0]} {value[1]} {value[2]}", name="reflectance")

# constant environment light
emitter = ET.SubElement(scene, "emitter", type="constant", id="World")
ET.SubElement(emitter, "rgb", value="1.000000 1.000000 1.000000", name="radiance")

# camera
sensor = ET.SubElement(scene, "sensor", type="perspective", id="Camera")
ET.SubElement(sensor, "string", name="fov_axis", value="x")
ET.SubElement(sensor, "float", name="fov", value=str(camera_settings["fov"]))
ET.SubElement(sensor, "float", name="principal_point_offset_x", value="0.000000")
ET.SubElement(sensor, "float", name="principal_point_offset_y", value="-0.000000")
ET.SubElement(sensor, "float", name="near_clip", value="0.100000")
ET.SubElement(sensor, "float", name="far_clip", value="10000.000000")

sionna_transform = ET.SubElement(sensor, "transform", name="to_world")
ET.SubElement(sionna_transform, "rotate", x="1", angle=str(camera_settings["rotation"][0]))
ET.SubElement(sionna_transform, "rotate", y="1", angle=str(camera_settings["rotation"][1]))
ET.SubElement(sionna_transform, "rotate", z="1", angle=str(camera_settings["rotation"][2]))

camera_position = np.array([0, 0, 100])  # Adjust camera height
ET.SubElement(sionna_transform, "translate", value=" ".join(map(str, camera_position)))

sampler = ET.SubElement(sensor, "sampler", type="independent")
ET.SubElement(sampler, "integer", name="sample_count", value="$spp")

film = ET.SubElement(sensor, "film", type="hdrfilm")
ET.SubElement(film, "integer", name="width", value="$resx")
ET.SubElement(film, "integer", name="height", value="$resy")



# 4. ground

#4326 bound
ground_polygon_4326 = shapely.geometry.Polygon(points)
ground_polygon_4326_bbox = ground_polygon_4326.bounds
print(ground_polygon_4326)


#utm bound
points_utm = [to_utm.transform(lon, lat) for lon, lat in points] 
ground_polygon_utm = shapely.geometry.Polygon(points_utm)
ground_polygon_utm_bbox = ground_polygon_utm.bounds

#envelop
ground_polygon_utm_envelop = ground_polygon_utm.envelope
center_x = ground_polygon_utm_envelop.centroid.x
center_y = ground_polygon_utm_envelop.centroid.y
print(ground_polygon_utm.exterior)

#vertex and edges
envelop_xy = unique_coords(reorder_localize_coords(ground_polygon_utm.exterior, center_x, center_y))
holes_xy = []

def edge_idx(len):
    arr = np.append(np.arange(len), 0)
    return np.stack([arr[:-1], arr[1:]], axis = 1)

v = []
e = []
nv = 0
for loop in (envelop_xy, *holes_xy):
    print(f"Loop: {loop}")
    v.append(loop)
    e.append(nv + edge_idx(len(loop)))
    nv += len(loop)

verts = np.concatenate(v)
edges = np.concatenate(e)

holes = np.array([np.mean(h, axis=0) for h in holes_xy])

#triangulate
result = triangulate(dict(vertices = verts, segments = edges), opts = "p")

#points and idx
vertices, indices = result["vertices"], result["triangles"]
nv, nf = len(vertices), len(indices)
vertices = np.concatenate([vertices, np.zeros((nv, 1))], axis=1)
print(f"points from triangulate: {vertices}" )

#write mesh
mesh = o3d.t.geometry.TriangleMesh()
mesh.vertex.positions = o3d.core.Tensor(vertices)
mesh.triangle.indices = o3d.core.Tensor(indices)

write_ply_ascii = False
mesh.scale(ground_scale, mesh.get_center())
o3d.t.io.write_triangle_mesh(mesh_dir / "ground.ply", mesh, write_ascii = write_ply_ascii)

#write xml
material_type = "mat-itu_medium_dry_ground"
sionna_shape = ET.SubElement(scene, "shape", type="ply", id=f"mesh-ground")
ET.SubElement(sionna_shape, "string", name="filename", value=f"test_mesh/ground.ply")
bsdf_ref = ET.SubElement(sionna_shape, "ref", id=material_type, name="bsdf")
ET.SubElement(sionna_shape, "boolean", name="face_normals", value="true")



# 5. select buildings

west = ground_polygon_4326_bbox[0]  # minx
south = ground_polygon_4326_bbox[1]  # miny
east = ground_polygon_4326_bbox[2]  # maxx
north = ground_polygon_4326_bbox[3]  # maxy

# width height
width = math.ceil(ground_polygon_utm_bbox[2] - ground_polygon_utm_bbox[0])
height = math.ceil(ground_polygon_utm_bbox[3] - ground_polygon_utm_bbox[1])
print(f"Estimated ground polygon size: width={width}m, height={height}m")
if width > 5000 or height > 5000:
    print(f"Too large!")
    exit(-1)

# Select buildings from the bundled OSM extract by default. Set
# RAYLOC_USE_OVERPASS=1 to refresh them from the public Overpass API.
osm_file = project_root / "assets" / "uw.osm"
if os.environ.get("RAYLOC_USE_OVERPASS") == "1":
    buildings = ox.features.features_from_bbox(
        bbox=ground_polygon_4326_bbox,
        tags={"building": True},
    )
else:
    # JOSM exports can retain deleted nodes without coordinates. OSMnx expects
    # every node to have latitude and longitude, so filter those records in a
    # temporary copy while leaving the bundled source extract unchanged.
    osm_tree = ET.parse(osm_file)
    osm_root = osm_tree.getroot()
    invalid_elements = [
        element for element in osm_root
        if element.get("visible") == "false"
        or element.get("action") == "delete"
        or (
            element.tag == "node"
            and (element.get("lat") is None or element.get("lon") is None)
        )
        or (element.tag == "way" and not element.findall("nd"))
    ]
    for element in invalid_elements:
        osm_root.remove(element)

    valid_way_ids = {way.get("id") for way in osm_root.findall("way")}
    invalid_relations = [
        relation for relation in osm_root.findall("relation")
        if any(
            member.get("type") == "way"
            and member.get("ref") not in valid_way_ids
            for member in relation.findall("member")
        )
    ]
    for relation in invalid_relations:
        osm_root.remove(relation)

    with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as temp_osm:
        sanitized_osm_file = Path(temp_osm.name)
    try:
        osm_tree.write(sanitized_osm_file, encoding="utf-8", xml_declaration=True)
        buildings = ox.features.features_from_xml(
            sanitized_osm_file,
            polygon=ground_polygon_4326,
            tags={"building": True},
        )
    finally:
        sanitized_osm_file.unlink(missing_ok=True)
buildings = buildings.to_crs(EPSG_code)

#remove intersected buildings
filtered_buildings = buildings[buildings.intersects(ground_polygon_utm)]
buildings_list = filtered_buildings.to_dict('records')

# 6. image
map_image = Image.new('L', (width, height), 0)



# 7.



# 8. write building
logger = logging.getLogger(__name__)

for idx, building in tqdm(enumerate(buildings_list), total=len(buildings_list), desc="Parsing buildings"):
    building_polygon = shape(building['geometry'])
    #print(building['geometry'])

    if (building_polygon.geom_type != 'Polygon'):
        logger.debug(f"building_polygon.geom_type: {building_polygon.geom_type}")
        continue

    building_height = random_building_height(building, building_polygon) #get height

    envelop_xy = unique_coords(reorder_localize_coords(building_polygon.exterior, center_x, center_y))

    holes_xy = []
    if (len(list(building_polygon.interiors)) != 0):
        for hole in list(building_polygon.interiors):
            valid_hole = reorder_localize_coords(hole, center_x, center_y)
            holes_xy.append(valid_hole)

    def edge_idx(len):
        arr = np.append(np.arange(len), 0)
        return np.stack([arr[:-1], arr[1:]], axis = 1)

    v = []
    e = []
    nv = 0
    for loop in (envelop_xy, *holes_xy):
        #print(f"Loop: {loop}")
        v.append(loop)
        e.append(nv + edge_idx(len(loop)))
        nv += len(loop)

    verts = np.concatenate(v)
    edges = np.concatenate(e)

    holes = np.array([np.mean(h, axis=0) for h in holes_xy])

    #triangulate
    if (len(holes) != 0):
        result = triangulate(dict(vertices = verts, segments = edges, holes = holes), opts = "p")
    else:
        result = triangulate(dict(vertices = verts, segments = edges), opts = "p")

    # print(verts)
    # print(edges)

    #points and idx
    vertices, indices = result["vertices"], result["triangles"]
    nv, nf = len(vertices), len(indices)
    vertices = np.concatenate([vertices, np.zeros((nv, 1))], axis=1)
    #print(f"points from triangulate: {vertices}" )

    #write mesh
    mesh = o3d.t.geometry.TriangleMesh()
    mesh.vertex.positions = o3d.core.Tensor(vertices)
    mesh.triangle.indices = o3d.core.Tensor(indices)

    write_ply_ascii = False
    wedge = mesh.extrude_linear([0, 0, building_height]) #apply height
    o3d.t.io.write_triangle_mesh(mesh_dir / f"building_{idx}.ply", wedge, write_ascii = write_ply_ascii)

    #write xml
    material_type = "mat-itu_concrete"
    sionna_shape = ET.SubElement(scene, "shape", type="ply", id=f"mesh-building_{idx}")
    ET.SubElement(sionna_shape, "string", name="filename", value=f"test_mesh/building_{idx}.ply")
    bsdf_ref = ET.SubElement(sionna_shape, "ref", id=material_type, name="bsdf")
    ET.SubElement(sionna_shape, "boolean", name="face_normals", value="true")
    
    #image
    local_exterior = reorder_localize_coords(building_polygon.exterior, ground_polygon_utm_envelop.bounds[0], ground_polygon_utm_envelop.bounds[3]) #for drawing, different offset
    ImageDraw.Draw(map_image).polygon([(x, -y) for x, y in list(local_exterior)], #plot -y because of display
                                            outline=int(building_height), fill=int(building_height))
    # print(ground_polygon_utm_envelop)
    # print(ground_polygon_utm_envelop.bounds)
    # print(center_x, center_y)
    # print(ground_polygon_utm_envelop.bounds[0], ground_polygon_utm_envelop.bounds[3])
    
    #break
    
xml_string = ET.tostring(scene, encoding="utf-8")
xml_pretty = minidom.parseString(xml_string).toprettyxml(indent="    ") 

with open(asset_dir / "scene.xml", "w", encoding="utf-8") as xml_file:
    xml_file.write(xml_pretty)

np.save(asset_dir / "2D_Building_Height_Map.npy", np.array(map_image))

plt.imshow(map_image)
plt.axis("off")
plt.savefig(asset_dir / "building_height_map.png", dpi=200, bbox_inches="tight")
plt.close()