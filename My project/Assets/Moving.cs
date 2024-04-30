using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Moving : MonoBehaviour
{
    public float height, speed, HForce = 2f,Rot = 0f;
    public Rigidbody RB;

    // Start is called before the first frame update
    void Start()
    {
        speed = 0f;
        height = 0f;
        RB = gameObject.GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    private void FixedUpdate()
    {
        if (gameObject.transform.position.y != height)
        {
            RB.AddForce(new Vector3(0, (height-gameObject.transform.position.y)*HForce * Time.fixedDeltaTime, 0));
            if(RB.velocity.y > Mathf.Abs(height - gameObject.transform.position.y))
            {
                RB.velocity = new Vector3(RB.velocity.x, Mathf.Abs(height - gameObject.transform.position.y), RB.velocity.z);
            }
        }
        if (Input.GetKey(KeyCode.UpArrow))
        {
            RB.AddRelativeForce(new Vector3(0,0, 1 * speed * Time.fixedDeltaTime),ForceMode.VelocityChange);
        }
        if (Input.GetKey(KeyCode.DownArrow))
        {
            RB.AddRelativeForce(new Vector3(0, 0, -1 * speed * Time.fixedDeltaTime), ForceMode.VelocityChange);
        }
        RB.rotation = Quaternion.Euler(0,Rot,0);
        if (RB.velocity.magnitude > 0.01) { RB.velocity -= RB.velocity.normalized * 0.5f * Time.fixedDeltaTime; } else
        {
            RB.velocity = Vector3.zero;
        }
    
    }
    void Update()
    {

        if (Input.GetKey(KeyCode.LeftArrow))
        {
            Rot += -3 * Time.deltaTime;

        }
        if (Input.GetKey(KeyCode.RightArrow))
        {
            Rot += 3 * Time.deltaTime;
        }
    }
}
